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
#include <COLLADAFWFileInfo.h>

#include <Math/COLLADABUMathMatrix4.h>
#include <Math/COLLADABUMathMatrix3.h>
#include <Math/COLLADABUMathQuaternion.h>
#include <Math/COLLADABUMathVector3.h>
#include <COLLADABUURI.h>

#include <Logging/Logger.h>
#include <string.h>

#define IN(s)                                \
    string oldspace = space;                 \
    space += "  ";                           \
    logger.info << space << s << logger.end;

#define OUT()                                \
    space = oldspace;

#define IN(s) 
#define OUT() 

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
  , upIndex(1)
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
    resource_dir = File::Parent(this->file);

    // logger.info << "Load Resource: " << file << logger.end;

    COLLADASaxFWL::Loader* loader = new COLLADASaxFWL::Loader();
    COLLADAFW::Root* root = new COLLADAFW::Root(loader, this);
    
    // Load scene graph 
    if (!root->loadDocument(file)) {
        Warning("Resource not loaded.");
        return;
    }
    if (!visualScene) {
        Warning("No visual scene found.");
        return;
    }
    this->root = new SceneNode();
    NodePointerArray& nodes = visualScene->getRootNodes();
    for (unsigned int i = 0; i < nodes.getCount(); i++) {
        this->root->AddNode(ReadNode(nodes[i]));
    }
    // logger.info << "Resource loaded" << logger.end;    
}

/**
 * Unload the resource.
 * Resets the root node. Does not delete the scene graph.
 */
void ColladaResource::Unload() {
    // delete (hopefully) all the intermediate structures ...
    for (map<COLLADAFW::UniqueId, GeoPrimitives*>::iterator i = geometries.begin();
         i != geometries.end(); 
         ++i) {
        for (GeoPrimitives::iterator j = (*i).second->begin();
             j != (*i).second->end();
             ++j) {
            delete *j;
        }
        delete (*i).second;
    }
    effects.clear();
    images.clear();
    geometries.clear();
    materials.clear();
    nodes.clear();
    if (visualScene) {
        delete visualScene;
        visualScene = NULL;
    }
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
    if (!m) Error("Invalid material id");
    MaterialPtr mp = effects[m->getInstantiatedEffect()];
    if (!mp) Error("Invalid effect id");
    return mp;
}

ISceneNode* ColladaResource::ReadTransformation(Transformation* t) {
    IN("+ReadTransformation");
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
        Warning("Unsupported transformation type: LOOKAT");
        break;
    case Transformation::SKEW:
        Warning("Unsupported transformation type: SKEW");
        break;
    default:
        Warning("Unsupported transformation type.");
    };
    OUT();
    return tn;
}
    
ISceneNode* ColladaResource::ReadInstanceGeometry(COLLADAFW::InstanceGeometry* ig) {
    IN("+ReadInstanceGeometry");
    GeoPrimitives* gps = LookupGeometry(ig->getInstanciatedObjectId());
    map<MaterialId, UniqueId> bindings = ExtractMaterialBindingMap(ig->getMaterialBindings());
    ISceneNode* sn = CreateGeometry(gps, bindings);
    OUT();
    return sn;
}

ColladaResource::GeoPrimitives* ColladaResource::ReadGeometry(const COLLADAFW::Geometry* g) {
    IN("+ReadGeometry");
    if (g->getType() != COLLADAFW::Geometry::GEO_TYPE_MESH) {
        Warning("Unsupported geometry type.");
        return NULL;
    }

    GeoPrimitives* gps = new GeoPrimitives();
    Mesh* mesh = (Mesh*)g;

    const int stride = 3;  // position and normal stride.
    float* posArray;
    bool delPos = ExtractFloatArray(mesh->getPositions(), &posArray);

    MeshVertexData& norm = mesh->getNormals();
    float* normArray;
    bool delNorm = ExtractFloatArray(norm, &normArray);

    MeshVertexData& uv = mesh->getUVCoords();
    float* uvArray;
    bool delUV = ExtractFloatArray(uv, &uvArray);
    MeshVertexData::InputInfos* uvInfo = ExtractInputInfos(uv);
    unsigned int  uvStride = uvInfo ? uvInfo->mStride : 0;
    // logger.info << space << "uvinputinfocount: " << uv.getNumInputInfos() << logger.end;
    // logger.info << space << "uvstride: " << uvInfo.mStride << " name: " << uvInfo.mName << logger.end;
    
    MeshVertexData& col = mesh->getColors();
    float* colArray;
    bool delCol = ExtractFloatArray(col, &colArray);
    
    MeshPrimitiveArray& prims = mesh->getMeshPrimitives();
    for (unsigned int i = 0; i < prims.getCount(); i++) {
        MeshPrimitive*    prim = prims[i];
        MaterialId mId = prim->getMaterialId(); 
        FaceSet* fs = new FaceSet();
        gps->push_back(new GeoPrimitive(mId, fs));

        unsigned int     count = prim->getFaceCount();
        UIntValuesArray&  posI = prim->getPositionIndices();
        UIntValuesArray& normI = prim->getNormalIndices();
        unsigned int*     colI = NULL;
        unsigned int*      uvI = NULL;
        // unsigned int indexStride = 0;        
        // if the primitive has some texture coordinate lists we
        // simply choose the first one.
        if (prim->hasUVCoordIndices()) {
            IndexList* uvIL = prim->getUVCoordIndices(0);
            uvI = uvIL->getIndices().getData();
            // logger.info << space << "uvindexstride: " << uvIL->getStride() <<  ". uvindexsetindex: " << uvIL->getSetIndex() << ". uvindexinitialindex: " << uvIL->getInitialIndex()  << logger.end;
            //logger.info << "posISize: " << 
            // indexStride = uvIL->getStride();
        }
        if (prim->hasColorIndices()) {
            IndexList* colIL = prim->getColorIndices(0);
            colI = colIL->getIndices().getData();
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
                    Vector<4,float> cols[3];
                    // for each vertex.
                    for (int k = 2; k >= 0; --k, ++index){
                        // cols[k][3] = 1.0;
                        // for each vertex component.
                        for (int l = 0; l < stride; l++) {
                            verts[k][l] = posArray[l+posI[index]*stride];
                            if (mesh->hasNormals()) 
                                norms[k][l] = normArray[l+normI[index]*stride];
                            if (colI)
                                cols[k][l] = colArray[l+colI[index]*stride];
                        }
                        if (uvI) {
                            for (unsigned int l = 0; l < 2; l++) {
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
                        if (colI) {
                            face->colr[0] = cols[0];
                            face->colr[1] = cols[1];
                            face->colr[2] = cols[2];
                            // face->colr[3] = cols[3];
                        }
                        //logger.info << "cols: " << face->colr[0] << logger.end;
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
            Warning("Ignoring unsupported primitive type.");
        };
    } 

    if (delPos) delete[] posArray;
    if (delNorm) delete[] normArray;
    if (delUV) delete[] uvArray;
    if (delCol) delete[] uvArray;
    OUT();
    return gps;
}

ISceneNode* ColladaResource::ReadNode(Node* node) {
    IN("+ReadNode");
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
    OUT();
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
            for (unsigned int k = 0; k < 3; ++k) {
                float tmp = f->vert[k][1];
                f->vert[k][1] = f->vert[k][upIndex];
                f->vert[k][upIndex] = tmp;
                tmp = f->norm[k][1];
                f->norm[k][1] = f->norm[k][upIndex];
                f->norm[k][upIndex] = tmp;
            }
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

bool ColladaResource::ExtractFloatAttribute(FloatOrParam& fop, float& dest) {
    if (fop.getType() == FloatOrParam::FLOAT && fop.getFloatValue() != -1) {
        dest = fop.getFloatValue();
        return true;
    }
    return false;
}

bool ColladaResource::ExtractFloatArray(MeshVertexData& d, float** dest) {
    *dest = NULL;
    if (d.empty()) 
        return false;
    unsigned int size = d.getValuesCount();
    switch (d.getType()) {
    case FloatOrDoubleArray::DATA_TYPE_FLOAT:
        *dest = d.getFloatValues()->getData();
        // memcpy(out, d.getFloatValues()->getData(), size * sizeof(float));
        break;
    case FloatOrDoubleArray::DATA_TYPE_DOUBLE: 
        {
            double* dv = d.getDoubleValues()->getData();
            *dest = new float[size];
            for (unsigned int i = 0; i < size; ++i) {
                *dest[i] = dv[i];
            };
            return true;
        }
    default:
        Warning("Unknown data array format.");
    }; 
    return false;
}

MeshVertexData::InputInfos* ColladaResource::ExtractInputInfos(MeshVertexData& d) {
    MeshVertexData::InputInfos* out;
    if (d.getNumInputInfos() == 0)
        return NULL;
    out = d.getInputInfosArray()[0];
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
    switch (asset->getUpAxisType()) {
    case FileInfo::X_UP:
        upIndex = 0;
        break;
    case FileInfo::Y_UP:
        upIndex = 1;
        break;
    case FileInfo::Z_UP:
        upIndex = 2;
        break;
    case FileInfo::NONE:
    default:
        Warning("No up-axis defined. Assuming y is up");
    };
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

string PrintWrap(Sampler::WrapMode m) {
    switch (m) {
    case Sampler::WRAP_MODE_NONE: 
        return "WM_NONE";
    case Sampler::WRAP_MODE_WRAP: 
        return "WM_WRAP";
    case Sampler::WRAP_MODE_MIRROR: 
        return "WM_MIRROR";
    case Sampler::WRAP_MODE_CLAMP: 
        return "WM_CLAMP";
    case Sampler::WRAP_MODE_BORDER: 
        return "WM_BORDER";
    default:
        return "WM_UNSPECIFIED";
    }
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
    ExtractFloatAttribute(ce->getShininess(), m->shininess);
    
    if (ce->getSamplerPointerArray().getCount() > 0) {
        logger.info << "sampler is present" << logger.end;
        Sampler* s = ce->getSamplerPointerArray()[0];
        if (s->getSamplerType() != Sampler::SAMPLER_TYPE_2D) {
            Warning("Unsupported texture sampling type");
        }
        else {
            //logger.info << "adding texture to material" << logger.end;
            //logger.info << "wrapS: " << PrintWrap(s->getWrapS()) << " wrapT: " << PrintWrap(s->getWrapT()) << logger.end;
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
    IN("writeImage");
    if (image->getSourceType() != Image::SOURCE_TYPE_URI) 
        Error("Unsupported image source type.");
    const URI& uri = image->getImageURI();
    //logger.info << "orguri: " << uri.originalStr() << " resuri: " << uri.getURIString() << logger.end;
    ITextureResourcePtr texr = ResourceManager<ITextureResource>::Create(resource_dir + uri.getURIString());
    images[image->getUniqueId()] = texr;
    OUT();
    return true;
}

    /** Writes the light.
		@return True on succeeded, false otherwise.*/
bool ColladaResource::writeLight( const COLLADAFW::Light* light ) {
    IN("writeLight");
    OUT();
    return true;
}

    /** Writes the animation.
		@return True on succeeded, false otherwise.*/
bool ColladaResource::writeAnimation( const COLLADAFW::Animation* animation ) {
    IN("writeAnimation");

    OUT();
    return true;
}

    /** Writes the animation.
		@return True on succeeded, false otherwise.*/
bool ColladaResource::writeAnimationList( const COLLADAFW::AnimationList* animationList ) {
    IN("writeAnimationList");
    OUT();
    return true;
}

    /** Writes the skin controller data.
		@return True on succeeded, false otherwise.*/
bool ColladaResource::writeSkinControllerData( const COLLADAFW::SkinControllerData* skinControllerData ) {
    IN("writeSkinControllerData");
    
    OUT();
    return true;
}

    /** Writes the controller.
		@return True on succeeded, false otherwise.*/
bool ColladaResource::writeController( const COLLADAFW::Controller* Controller ) {
    IN("writeController");

    OUT();
    return true;
}

    /** When this method is called, the writer must write the formulas. All the formulas of the entire
		COLLADA file are contained in @a formulas.
		@return The writer should return true, if writing succeeded, false otherwise.*/
bool ColladaResource::writeFormulas( const COLLADAFW::Formulas* formulas ) {
    IN("writeFormulas");

    OUT();
    return true;
}

    /** When this method is called, the writer must write the kinematics scene. 
		@return The writer should return true, if writing succeeded, false otherwise.*/
bool ColladaResource::writeKinematicsScene( const COLLADAFW::KinematicsScene* kinematicsScene ) {
    IN("writeKinematics");

    OUT();
    return true;
}

} // NS Resources
} // NS OpenEngine
