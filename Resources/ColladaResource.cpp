// Collada model resource.
// -------------------------------------------------------------------
// Copyright (C) 2007 OpenEngine.dk (See AUTHORS) 
// 
// 
// This program is free software; It is covered by the GNU General 
// Public License version 2 or any later version. 
// See the GNU General Public License for more details (see LICENSE). 
//--------------------------------------------------------------------

#include <Resources/ColladaResource.h>

#include <Scene/SceneNode.h>
#include <Geometry/Material.h>
#include <Geometry/FaceSet.h>
#include <Scene/ISceneNode.h>
#include <Scene/SceneNode.h>
#include <Scene/TransformationNode.h>
#include <Scene/GeometryNode.h>
#include <Math/Matrix.h>
#include <Math/Quaternion.h>
#include <Resources/Exceptions.h>
#include <Resources/ResourceManager.h>
#include <Resources/File.h>

// OpenCollada stuff
#include <COLLADASaxFWLLoader.h>
#include <COLLADAFWRoot.h>
#include <COLLADAFWVisualScene.h>
#include <COLLADAFWEffectCommon.h>
#include <COLLADAFWEffect.h>
#include <COLLADAFWColorOrTexture.h>
#include <COLLADAFWFloatOrParam.h>
#include <COLLADAFWColor.h>
#include <COLLADAFWMaterial.h>
#include <COLLADAFWGeometry.h>
#include <COLLADAFWMesh.h>
#include <COLLADAFWMeshPrimitive.h>
#include <COLLADAFWTransformation.h>
#include <COLLADAFWMatrix.h>
#include <COLLADAFWTranslate.h>
#include <COLLADAFWRotate.h>
#include <COLLADAFWScale.h>
#include <COLLADAFWLibraryNodes.h>
#include <COLLADAFWImage.h>

#include <Math/COLLADABUMathMatrix4.h>
#include <Math/COLLADABUMathMatrix3.h>
#include <Math/COLLADABUMathQuaternion.h>
#include <Math/COLLADABUMathVector3.h>
#include <COLLADABUURI.h>

#include <Logging/Logger.h>

#include <string.h>


namespace OpenEngine {
namespace Resources {

using namespace Scene;
using namespace Geometry;

using namespace COLLADAFW;
using namespace COLLADABU::Math;
using namespace COLLADABU;

/**
 * Get the file extension for Collada files.
 */
ColladaPlugin::ColladaPlugin() {
    this->AddExtension("dae");
}

/**
 * Create a Collada resource.
 */
IModelResourcePtr ColladaPlugin::CreateResource(string file) {
    return IModelResourcePtr(new ColladaResource(file));
}


/**
 * Resource constructor.
 */
ColladaResource::ColladaResource(string file)
  : file(file)
  , root(NULL)
  , visualScene(NULL) 
{}

/**
 * Resource destructor.
 */
ColladaResource::~ColladaResource() {
    Unload();
}

/**
 * Load a Collada dae scene file.
 *
 * This method parses the file given to the constructor and populates a
 * scene graph with the data from the file that can be retrieved with
 * GetSceneNode().
 *
 * @see Scene::ISceneNode
 */
void ColladaResource::Load() {
    if (root) return;

    // we reset the resource path temporary to create the texture resource
    string resource_dir = File::Parent(this->file);

	if (! DirectoryManager::IsInPath(resource_dir)) {
        DirectoryManager::AppendPath(resource_dir);
    }

    logger.info << "Load Resource: " << file << logger.end;

    COLLADASaxFWL::Loader* loader = new COLLADASaxFWL::Loader();
    COLLADAFW::Root* root = new COLLADAFW::Root(loader, this);
    
    // Load scene graph 
    if (!root->loadDocument(file)) {
        logger.info << "Resource not loaded." << logger.end;
        return;
    }
    if (!visualScene) {
        logger.warning << "No visual scene found." << logger.end;
        return;
    }
    this->root = new SceneNode();
    NodePointerArray& nodes = visualScene->getRootNodes();
    for (unsigned int i = 0; i < nodes.getCount(); i++) {
        this->root->AddNode(ReadNode(nodes[i]));
    }
    logger.info << "Resource loaded" << logger.end;    
}

/**
 * Unload the resource.
 * Resets the root node. Does not delete the scene graph.
 */
void ColladaResource::Unload() {
    effects.clear();
    images.clear();
    geometries.clear();
    materials.clear();
    nodes.clear();
    // if (visualScene) 
    //     delete visualScene;
    root = NULL;
}

/**
 * Get the scene graph for the loaded collada data.
 *
 * @return Root node of the scene graph
 */
ISceneNode* ColladaResource::GetSceneNode() {
    return root;
}

void ColladaResource::Error(string msg) {
    logger.error << "OpenCollada: " << msg << logger.end;
    throw new ResourceException("OpenCollada: " + msg);
}

void ColladaResource::Warning(string msg) {
    logger.warning << "OpenCollada: " <<  msg << logger.end;
}

Node* ColladaResource::LookupNode(UniqueId id) {
    Node* n = nodes[id];
    if (!n) Error("Invalid node id");
    return n;
}

ColladaResource::GeoPrimitives* ColladaResource::LookupGeometry(UniqueId id) {
    GeoPrimitives* gps = geometries[id];
    if (!gps) Error("Invalid geometry id");
    return gps;
}

ITextureResourcePtr ColladaResource::LookupImage(UniqueId id) {
    ITextureResourcePtr texr = images[id];
    if (!texr) Error("Invalid image id");
    return texr;
}

MaterialPtr ColladaResource::LookupMaterial(UniqueId id) {
    const COLLADAFW::Material* m = materials[id];
    logger.info << "lookup material: " << m->getName() << logger.end;
    if (!m) Error("Invalid material id");
    MaterialPtr mp = effects[m->getInstantiatedEffect()];
    if (!mp) Error("Invalid effect id");
    return mp;
}

ISceneNode* ColladaResource::ReadTransformation(Transformation* t) {
    string oldspace = space;
    space += "  ";
    logger.info << space << "+ReadTransformation" << logger.end;
    TransformationNode* tn = new TransformationNode(); 
    switch (t->getTransformationType()) {
    case Transformation::MATRIX: 
        {
            Matrix4      _m = ((COLLADAFW::Matrix*)t)->getMatrix();
            Vector3    _pos = _m.getTrans();
            Vector3    _scl = _m.getScale();
            Matrix3    _rot;
            _m.extract3x3Matrix(_rot);
            Matrix<3,3,float> m(_m[0][0], _m[0][1], _m[0][2],
                                _m[1][0], _m[1][1], _m[1][2],
                                _m[2][0], _m[2][1], _m[2][2]);
            tn->SetRotation(Quaternion<float>(m));
            tn->SetPosition(Vector<3,float>(_pos[0], _pos[1], _pos[2]));
            tn->SetScale(Vector<3,float>(_scl[0], _scl[1], _scl[2]));
            break;
        }
    case Transformation::TRANSLATE:
        {
            Vector3& _pos = ((Translate*)t)->getTranslation();
            tn->SetPosition(Vector<3,float>(_pos[0], _pos[1], _pos[2]));
        }        
        break;
    case Transformation::ROTATE:
        {
            Rotate* _rot = ((Rotate*)t);
            float     _a = _rot->getRotationAngle();
            Vector3&  _v = _rot->getRotationAxis();
            tn->SetRotation(Quaternion<float>(_a, Vector<3,float>(_v[0], _v[1], _v[2])));
        }
        break;
    case Transformation::SCALE:
        {
            Vector3&  _scl = ((Scale*)t)->getScale();
            tn->SetScale(Vector<3,float>(_scl[0], _scl[1], _scl[2]));
            break;
        }
    case Transformation::LOOKAT:
        logger.warning << "Unsupported transformation type: LOOKAT" << logger.end;
        break;
    case Transformation::SKEW:
        logger.warning << "Unsupported transformation type: SKEW" << logger.end;
        break;
    default:
        logger.warning << "Unsupported transformation type." << logger.end;
    };
    space = oldspace;
    return tn;
}
    
ISceneNode* ColladaResource::ReadInstanceGeometry(COLLADAFW::InstanceGeometry* ig) {
    string oldspace = space;
    space += "  ";
    logger.info << space << "+ReadInstanceGeometry" << logger.end;
    GeoPrimitives* gps = LookupGeometry(ig->getInstanciatedObjectId());
    map<MaterialId, UniqueId> bindings = ExtractMaterialBindingMap(ig->getMaterialBindings());
    ISceneNode* sn = CreateGeometry(gps, bindings);
    space = oldspace;
    return sn;
}

ColladaResource::GeoPrimitives* ColladaResource::ReadGeometry(const COLLADAFW::Geometry* g) {
    string oldspace = space;
    space += "  ";
    logger.info << space << "+ReadGeometry" << logger.end;


    if (g->getType() != COLLADAFW::Geometry::GEO_TYPE_MESH) {
        logger.warning << "Unsupported geometry type." << logger.end;
        return NULL;
    }

    GeoPrimitives* gps = new GeoPrimitives();
    Mesh* mesh = (Mesh*)g;

    const int stride = 3;  // position and normal stride.
    float* posArray = ExtractFloatArray(mesh->getPositions());

    MeshVertexData& norm = mesh->getNormals();
    float* normArray = ExtractFloatArray(norm);

    MeshVertexData& uv = mesh->getUVCoords();
    float* uvArray = ExtractFloatArray(uv);
    // MeshVertexData::InputInfos uvInfo = ExtractInputInfos(uv);
    // logger.info << space << "uvinputinfocount: " << uv.getNumInputInfos() << logger.end;
    // logger.info << space << "uvstride: " << uvInfo.mStride << " name: " << uvInfo.mName << logger.end;

    // MeshVertexData& col = mesh->getColors();
    // float* colArray = ExtractFloatArray(col);
    // MeshVertexData::InputInfos colInfo = ExtractInputInfos(col);

    MeshPrimitiveArray& prims = mesh->getMeshPrimitives();
    for (unsigned int i = 0; i < prims.getCount(); i++) {
        MeshPrimitive*    prim = prims[i];
        MaterialId mId = prim->getMaterialId(); 
        FaceSet* fs = new FaceSet();
        gps->push_back(new GeoPrimitive(mId, fs));

        unsigned int     count = prim->getFaceCount();
        UIntValuesArray&  posI = prim->getPositionIndices();
        UIntValuesArray& normI = prim->getNormalIndices();
        unsigned int*      uvI = NULL;
        unsigned int  uvStride = 2;
        // if the primitive has some texture coordinate lists we
        // simply choose the first one.
        if (prim->hasUVCoordIndices()) {
            IndexList* uvIL = prim->getUVCoordIndices(0);
            logger.info << space << "uvindexstride: " << uvIL->getStride() <<  ". uvindexsetindex: " << uvIL->getSetIndex() << ". uvindexinitialindex: " << uvIL->getInitialIndex()  << logger.end;
            uvI = uvIL->getIndices().getData();
            // uvStride = uvIL->getStride();
      }
        // unsigned int* colI = mp->getColorIndicesArray()[0]->getIndices();
        switch (prim->getPrimitiveType()) {
        case MeshPrimitive::TRIANGLES: 
            {
                unsigned int index = 0;
                // for each face.
                for (unsigned int j = 0; j < count; j++) {
                    // intermediate face data
                    Vector<3,float> verts[3];
                    Vector<3,float> norms[3];
                    Vector<2,float> texc[3];
                    // Vector<4,float> cols[3];
                    // for each vertex.
                    for (int k = 0; k < 3; k++, index++){
                        // for each vertex component.
                        for (int l = 0; l < stride; l++) {
                            verts[k][l] = posArray[l+posI[index]*stride];
                            if (mesh->hasNormals()) 
                                norms[k][l] = normArray[l+normI[index]*stride];
                            //cols[k][l] = colArray[l+colI[index]*stride];
                        }
                        if (uvI) {
                            for (int l = 0; l < uvStride; l++) {
                                texc[k][l] = uvArray[l+uvI[index]*uvStride];
                            }
                        }
                    }
                    try {
                        FacePtr face = FacePtr(new Face(verts[0], verts[1], verts[2]));
                        face->norm[0] = norms[0];
                        face->norm[1] = norms[1];
                        face->norm[2] = norms[2];
                        face->texc[0] = texc[0];
                        face->texc[1] = texc[1];
                        face->texc[2] = texc[2];
                        // logger.info << "texc: " << texc[0] << logger.end;
                        // face->mat = mat;
                        fs->Add(face);
                    }
                    catch (Exception e) {
                        Warning("Disregarding invalid face.");
                    }
                }
                break;
            }
        default:
            logger.warning << "Ignoring unsupported primitive type." << logger.end;
        };
    } 


    // FaceSet* fs = new FaceSet();    

    space = oldspace;
    return gps;
}

ISceneNode* ColladaResource::ReadNode(Node* node) {
    string oldspace = space;
    space += "  ";
    logger.info << space << "+ReadNode" << logger.end;
    ISceneNode* r = new SceneNode();
    ISceneNode* c = r;

    // Read transformations
    TransformationNode* tn = new TransformationNode();
    Matrix4   _m = node->getTransformationMatrix();
    Vector3 _pos = _m.getTrans();
    Vector3 _scl = _m.getScale();
    Matrix3 _rot;
    _m.extract3x3Matrix(_rot);
    Matrix<3,3,float> m(_m[0][0], _m[0][1], _m[0][2],
                        _m[1][0], _m[1][1], _m[1][2],
                        _m[2][0], _m[2][1], _m[2][2]);
    tn->SetRotation(Quaternion<float>(m));
    tn->SetPosition(Vector<3,float>(_pos[0], _pos[1], _pos[2]));
    tn->SetScale(Vector<3,float>(_scl[0], _scl[1], _scl[2]));
    c->AddNode(tn);
    c = tn;
    // for (unsigned int i = 0; i < node->getTransformations().getCount(); i++) {
    //     ISceneNode* n = ReadTransformation(node->getTransformations()[i]);
    //     c->AddNode(n);
    //     c = n;
    // }

    // Read instance geometry
    for (unsigned int i = 0; i < node->getInstanceGeometries().getCount(); i++) {
        ISceneNode* n = ReadInstanceGeometry(node->getInstanceGeometries()[i]);
        c->AddNode(n);
    }

    // Read sub-nodes
    for (unsigned int i = 0; i < node->getInstanceNodes().getCount(); i++) {
        ISceneNode* n = ReadNode(LookupNode(node->getInstanceNodes()[i]->getInstanciatedObjectId()));
        c->AddNode(n);
    }
    for (unsigned int i = 0; i < node->getChildNodes().getCount(); i++) {
        ISceneNode* n = ReadNode(node->getChildNodes()[i]);
        c->AddNode(n);
    }
    space = oldspace;
    return r;
}

ISceneNode* ColladaResource::CreateGeometry(GeoPrimitives* gps, map<MaterialId, UniqueId> bindings) {
    SceneNode* sn = new SceneNode();
    for (GeoPrimitives::iterator i = gps->begin(); i != gps->end(); ++i) {
        GeoPrimitive* gp = *i;
        MaterialPtr m = LookupMaterial(bindings[gp->mId]);
        FaceSet* fs = new FaceSet();
        for (FaceList::iterator j = gp->fs->begin(); j != gp->fs->end(); ++j) {
            FacePtr f = FacePtr(new Face(*(*j)));
            f->mat = m;
            fs->Add(f);
        }
        GeometryNode* gn = new GeometryNode(fs);
        sn->AddNode(gn);
    }
    return sn;
}

map<MaterialId, UniqueId> ColladaResource::ExtractMaterialBindingMap(MaterialBindingArray& mbs) {
    map<MaterialId, UniqueId> bindings;
    for (unsigned int i = 0; i < mbs.getCount(); i++) {
        MaterialBinding b = mbs[i];
        bindings[b.getMaterialId()] = b.getReferencedMaterial();
    }
    return bindings;
}

bool ColladaResource::ExtractColor(ColorOrTexture& cot, Vector<4,float>& dest) {
    if (cot.isColor()) {
        Color& c = cot.getColor();
        dest[0] = c.getRed();
        dest[1] = c.getGreen();
        dest[2] = c.getBlue();
        dest[3] = c.getAlpha();
        return true;
    }
    return false;
}

bool ColladaResource::ExtractFloat(FloatOrParam& fop, float& dest) {
    if (fop.getType() == FloatOrParam::FLOAT) {
        dest = fop.getFloatValue();
        return true;
    }
    return false;
}

float* ColladaResource::ExtractFloatArray(MeshVertexData& d) {
    if (d.empty())
        return NULL;
    unsigned int size = d.getValuesCount();
    float* out = new float[size];
    double* dv = NULL;
    switch (d.getType()) {
    case FloatOrDoubleArray::DATA_TYPE_FLOAT:
        memcpy(out, d.getFloatValues()->getData(), size * sizeof(float));
        break;
    case FloatOrDoubleArray::DATA_TYPE_DOUBLE:
        dv = d.getDoubleValues()->getData();
        for (unsigned int i = 0; i < size; ++i) {
            out[i] = dv[i];
        };
        break;
    default:
        logger.warning << "Unknown data array format." << logger.end;
        delete out;
        out = NULL;
    }; 
    return out;
}

MeshVertexData::InputInfos ColladaResource::ExtractInputInfos(MeshVertexData& d) {
    MeshVertexData::InputInfos out;
    out.mName = "";
    out.mStride = 1;
    out.mLength = 0;
    if (d.getNumInputInfos() == 0)
        return out;
    MeshVertexData::InputInfos* _out = d.getInputInfosArray()[0];
    out.mName   = _out->mName;
    out.mStride = _out->mStride;
    out.mLength = _out->mLength;
    return out;
}

/** Deletes the entire scene.
	@param errorMessage A message containing informations about the error that occurred.
*/
void ColladaResource::cancel(const string& errorMessage) {
    Error(errorMessage);
}

/** Prepare to receive data.*/
void ColladaResource::start() {
    logger.info << "start loading" << logger.end;
}

/** Remove all objects that don't have an object. Deletes unused visual scenes.*/
void ColladaResource::finish() {
    logger.info << "finish loading" << logger.end;
}

/** When this method is called, the writer must write the global document asset.
	@return The writer should return true, if writing succeeded, false otherwise.*/
bool ColladaResource::writeGlobalAsset ( const COLLADAFW::FileInfo* asset ) {
    logger.info << "Global Asset" << logger.end;
    return true;
}

/** Writes the entire visual scene.
    @return True on succeeded, false otherwise.*/
bool ColladaResource::writeVisualScene ( const COLLADAFW::VisualScene* visualScene ) {
    logger.info << "Visual Scene" << logger.end;
    this->visualScene = new COLLADAFW::VisualScene(*visualScene);
    return true;
}

/** Writes the scene.
	@return True on succeeded, false otherwise.*/
bool ColladaResource::writeScene ( const COLLADAFW::Scene* scene ) {
    logger.info << "Scene" << logger.end;
    return true;
}

/** Handles all nodes in the library nodes.
	@return True on succeeded, false otherwise.*/
bool ColladaResource::writeLibraryNodes( const COLLADAFW::LibraryNodes* libraryNodes ) {
    logger.info << "Library Nodes" << logger.end;
    for (unsigned int i = 0; i < libraryNodes->getNodes().getCount(); i++) {
        Node* n = libraryNodes->getNodes()[i];
        nodes[n->getUniqueId()] = n;
    }
    return true;
}

/** Writes the geometry.
	@return True on succeeded, false otherwise.*/
bool ColladaResource::writeGeometry ( const COLLADAFW::Geometry* geometry ) {
    logger.info << "Geometry" << logger.end;

    geometries[geometry->getUniqueId()] = ReadGeometry(geometry);
    
    


    return true;
}

/** Writes the material.
	@return True on succeeded, false otherwise.*/
bool ColladaResource::writeMaterial( const COLLADAFW::Material* material ) {
    logger.info << "Material" << logger.end;
    materials[material->getUniqueId()] = new COLLADAFW::Material(*material);
    return true;
}

    /** Writes the effect.
		@return True on succeeded, false otherwise.*/
bool ColladaResource::writeEffect( const COLLADAFW::Effect* effect ) {
    logger.info << "Effect" << logger.end;
    // effects[effect->getUniqueId()] = new COLLADAFW::Effect(*effect);

    CommonEffectPointerArray cepa = effect->getCommonEffects();
    if (cepa.getCount() == 0) {
        logger.warning << "Effect contains no common effects." << logger.end;
        return true;
    }
    EffectCommon* ce = cepa[0];
    MaterialPtr m = MaterialPtr(new Material);

    ExtractColor(ce->getAmbient(), m->ambient);
    ExtractColor(ce->getDiffuse(), m->diffuse);
    ExtractColor(ce->getSpecular(), m->specular);
    ExtractColor(ce->getEmission(), m->emission);
    ExtractFloat(ce->getShininess(), m->shininess);
    if (m->shininess < 0.0) m->shininess = 0.0;
    
    if (ce->getSamplerPointerArray().getCount() > 0) {
        logger.info << "sampler is present" << logger.end;
        Sampler* s = ce->getSamplerPointerArray()[0];
        if (s->getSamplerType() != Sampler::SAMPLER_TYPE_2D) {
            Warning("Unsupported texture sampling type");
        }
        else {
            logger.info << "adding texture to material" << logger.end;
            m->texr = LookupImage(s->getSourceImage());
        }
    }
    
    effects[effect->getUniqueId()] = m;

    // @todo: handle shader type
    return true;
}

    /** Writes the camera.
		@return True on succeeded, false otherwise.*/
bool ColladaResource::writeCamera( const COLLADAFW::Camera* camera ) {
    logger.info << "Camera" << logger.end;
    return true;
}

    /** Writes the image.
		@return True on succeeded, false otherwise.*/
bool ColladaResource::writeImage( const COLLADAFW::Image* image ) {
    logger.info << "Image" << logger.end;
    if (image->getSourceType() != Image::SOURCE_TYPE_URI) 
        Error("Unsupported image source type.");
    const URI& uri = image->getImageURI();
    logger.info << "orguri: " << uri.originalStr() << " resuri: " << uri.getURIString() << logger.end;
    ITextureResourcePtr texr = ResourceManager<ITextureResource>::Create(uri.getURIString());
    images[image->getUniqueId()] = texr;
    return true;
}

    /** Writes the light.
		@return True on succeeded, false otherwise.*/
bool ColladaResource::writeLight( const COLLADAFW::Light* light ) {
    logger.info << "Light" << logger.end;
    return true;
}

    /** Writes the animation.
		@return True on succeeded, false otherwise.*/
bool ColladaResource::writeAnimation( const COLLADAFW::Animation* animation ) {
    logger.info << "Animation" << logger.end;
    return true;
}

    /** Writes the animation.
		@return True on succeeded, false otherwise.*/
bool ColladaResource::writeAnimationList( const COLLADAFW::AnimationList* animationList ) {
    logger.info << "Animation List" << logger.end;
    return true;
}

    /** Writes the skin controller data.
		@return True on succeeded, false otherwise.*/
bool ColladaResource::writeSkinControllerData( const COLLADAFW::SkinControllerData* skinControllerData ) {
    logger.info << "Skin Controller" << logger.end;
    return true;
}

    /** Writes the controller.
		@return True on succeeded, false otherwise.*/
bool ColladaResource::writeController( const COLLADAFW::Controller* Controller ) {
    logger.info << "Controller" << logger.end;
    return true;
}

    /** When this method is called, the writer must write the formulas. All the formulas of the entire
		COLLADA file are contained in @a formulas.
		@return The writer should return true, if writing succeeded, false otherwise.*/
bool ColladaResource::writeFormulas( const COLLADAFW::Formulas* formulas ) {
    logger.info << "Formulas" << logger.end;
    return true;
}

    /** When this method is called, the writer must write the kinematics scene. 
		@return The writer should return true, if writing succeeded, false otherwise.*/
bool ColladaResource::writeKinematicsScene( const COLLADAFW::KinematicsScene* kinematicsScene ) {
    logger.info << "Kinematics" << logger.end;
    return true;
}

} // NS Resources
} // NS OpenEngine
