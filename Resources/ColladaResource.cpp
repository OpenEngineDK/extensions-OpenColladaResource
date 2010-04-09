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
#include <Math/Matrix.h>
#include <Math/Quaternion.h>
#include <Resources/Exceptions.h>
#include <Resources/ITexture2D.h>
#include <Resources/ResourceManager.h>
#include <Resources/File.h>

#include <Scene/MeshNode.h>
#include <Geometry/Mesh.h>
#include <Geometry/GeometrySet.h>
#include <Resources/DataBlock.h>


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
#include <COLLADAFWScene.h>
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

#define _IN(s)                                \
    string oldspace = space;                 \
    space += "  ";                           \
    logger.info << space << s << logger.end;

#define _OUT()                                \
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
  // , upIndex(1)
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

    total = 0;
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
        ReadNode(nodes[i], this->root);
    }
    logger.info << "Resource loaded, vertex count: " << total << logger.end;    
}

/**
 * Unload the resource.
 * Resets the root node. Does not delete the scene graph.
 */
void ColladaResource::Unload() {
    // delete (hopefully) all the intermediate structures ...
    // for (map<COLLADAFW::UniqueId, GeoPrimitives*>::iterator i = geometries.begin();
    //      i != geometries.end(); 
    //      ++i) {
    //     for (GeoPrimitives::iterator j = (*i).second->begin();
    //          j != (*i).second->end();
    //          ++j) {
    //         delete *j;
    //     }
    //     delete (*i).second;
    // }
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

// ISceneNode* ColladaResource::ReadTransformation(Transformation* t) {
//     IN("+ReadTransformation");
//     TransformationNode* tn = new TransformationNode(); 
//     switch (t->getTransformationType()) {
//     case Transformation::MATRIX: 
//         {
//             Matrix4      _m = ((COLLADAFW::Matrix*)t)->getMatrix();
//             Vector3    _pos = _m.getTrans();
//             Vector3    _scl = _m.getScale();
//             Matrix3    _rot;
//             _m.extract3x3Matrix(_rot);
//             Matrix<3,3,float> m(_m[0][0], _m[0][1], _m[0][2],
//                                 _m[1][0], _m[1][1], _m[1][2],
//                                 _m[2][0], _m[2][1], _m[2][2]);
//             tn->SetRotation(Quaternion<float>(m));
//             tn->SetPosition(Vector<3,float>(_pos[0], _pos[1], _pos[2]));
//             tn->SetScale(Vector<3,float>(_scl[0], _scl[1], _scl[2]));
//             break;
//         }
//     case Transformation::TRANSLATE:
//         {
//             Vector3& _pos = ((Translate*)t)->getTranslation();
//             tn->SetPosition(Vector<3,float>(_pos[0], _pos[1], _pos[2]));
//         }        
//         break;
//     case Transformation::ROTATE:
//         {
//             Rotate* _rot = ((Rotate*)t);
//             float     _a = _rot->getRotationAngle();
//             Vector3&  _v = _rot->getRotationAxis();
//             tn->SetRotation(Quaternion<float>(_a, Vector<3,float>(_v[0], _v[1], _v[2])));
//         }
//         break;
//     case Transformation::SCALE:
//         {
//             Vector3&  _scl = ((Scale*)t)->getScale();
//             tn->SetScale(Vector<3,float>(_scl[0], _scl[1], _scl[2]));
//             break;
//         }
//     case Transformation::LOOKAT:
//         Warning("Unsupported transformation type: LOOKAT");
//         break;
//     case Transformation::SKEW:
//         Warning("Unsupported transformation type: SKEW");
//         break;
//     default:
//         Warning("Unsupported transformation type.");
//     };
//     OUT();
//     return tn;
// }
    
void ColladaResource::ReadInstanceGeometry(COLLADAFW::InstanceGeometry* ig, ISceneNode* parent) {
    IN("+ReadInstanceGeometry");
    GeoPrimitives* gps = LookupGeometry(ig->getInstanciatedObjectId());
    map<MaterialId, UniqueId> bindings = ExtractMaterialBindingMap(ig->getMaterialBindings());
    ISceneNode* node = new SceneNode();
    for (GeoPrimitives::iterator i = gps->begin(); i != gps->end(); ++i) {
        GeoPrimitive gp = *i;
        MaterialPtr mat = LookupMaterial(bindings[gp.mId]);
        Mesh mesh = gp.prim;
        node->AddNode(new MeshNode(MeshPtr(new Mesh(mesh.GetIndices(), 
                                                    mesh.GetType(), 
                                                    mesh.GetGeometrySet(),
                                                    mat))));
    }
    parent->AddNode(node);
    OUT();
}

ColladaResource::GeoPrimitives* ColladaResource::ReadGeometry(const COLLADAFW::Geometry* g) {
    IN("+ReadGeometry");
    if (g->getType() != COLLADAFW::Geometry::GEO_TYPE_MESH) {
        Warning("Unsupported geometry type.");
        return NULL;
    }

    GeoPrimitives* gps = new GeoPrimitives();
    COLLADAFW::Mesh* mesh = (COLLADAFW::Mesh*)g;

    const int stride = 3;  // position and normal stride.
    float* posArray;
    bool delPos = ExtractFloatArray(mesh->getPositions(), &posArray);

    MeshVertexData& norm = mesh->getNormals();
    float* normArray;
    bool delNorm = ExtractFloatArray(norm, &normArray);

    MeshVertexData& uv = mesh->getUVCoords();
    float* uvArray;
    bool delUV = ExtractFloatArray(uv, &uvArray);
    
    MeshVertexData& col = mesh->getColors();
    float* colArray;
    bool delCol = ExtractFloatArray(col, &colArray);
    
    MeshPrimitiveArray& prims = mesh->getMeshPrimitives();
    for (unsigned int i = 0; i < prims.getCount(); i++) {
        MeshPrimitive* prim = prims[i];
        MaterialId mId = prim->getMaterialId(); 
        
        unsigned int count     = prim->getFaceCount();
        unsigned int vertCount = count * 3;
        total += vertCount;
        UIntValuesArray& posI  = prim->getPositionIndices();
        UIntValuesArray& normI = prim->getNormalIndices();
        unsigned int* uvI      = NULL;
        unsigned int uvOffset  = 0;
        unsigned int uvStride  = 0;
        unsigned int* colI     = NULL;
        unsigned int colOffset = 0;
        unsigned int colStride = 0;

        float* vsArr  = new float[vertCount * 3];
        float* nsArr  = new float[vertCount * 3];
        float* uvArr  = NULL;
        float* colArr = NULL;
        
        unsigned int*   isArr = new unsigned int[vertCount];
        Float3DataBlockPtr vs = Float3DataBlockPtr(new DataBlock<3,float>(vertCount, vsArr));
        Float3DataBlockPtr ns = Float3DataBlockPtr(new DataBlock<3,float>(vertCount, nsArr));
        IndicesPtr         is = IndicesPtr(new Indices(vertCount, isArr));

        IDataBlockList uvs; 
        Float2DataBlockPtr uv;
        Float3DataBlockPtr col;

        // if the primitive has multiple texture coordinates we simply
        // choose the first one.
        if (prim->hasUVCoordIndices()) {
            IndexList* uvIL = prim->getUVCoordIndices(0);
            uvI = uvIL->getIndices().getData();
            uvOffset = uvIL->getInitialIndex();
            uvStride = uvIL->getStride();
            uvArr = new float[vertCount * 2];
            uvs.push_back(Float2DataBlockPtr(new DataBlock<2,float>(vertCount, uvArr)));
        }
        if (prim->hasColorIndices()) {
            IndexList* colIL = prim->getColorIndices(0);
            colI = colIL->getIndices().getData();
            colOffset = colIL->getInitialIndex();
            colStride = colIL->getStride();
            colArr = new float[vertCount * 3];
            col = Float3DataBlockPtr(new DataBlock<3,float>(vertCount, colArr));
            // logger.info << "colorStride: " << colIL->getStride() << logger.end;
        }

        switch (prim->getPrimitiveType()) {
        case MeshPrimitive::TRIANGLES: 
            {
                gps->push_back(GeoPrimitive(mId, 
                                            Mesh(is,
                                                 TRIANGLES, 
                                                 GeometrySetPtr(new GeometrySet(vs, ns, uvs, col)), 
                                                 MaterialPtr())));
                unsigned int index = 0;
                // for each face.
                for (unsigned int j = 0; j < count; j++) {
                    // for each vertex.
                    for (int k = 0; k < 3; ++k, ++index){
                        isArr[index] = index;
                        for (int l = 0; l < stride; l++) {
                            vsArr[index*3+coord[l]] = posArray[l+posI[index]*stride];
                            if (mesh->hasNormals()) {
                                nsArr[index*3+coord[l]] = normArray[l+normI[index]*stride];
                            }
                            if (colI) {
                                colArr[index*3+l] = colArray[l+colI[index+colOffset]*colStride];
                                
                            }
                        }
                        if (uvI) {
                            for (unsigned int l = 0; l < 2; l++) {
                                uvArr[index*2+l] = uvArray[l+uvI[index]*uvStride];
                            }
                        }
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

void ColladaResource::ReadNode(Node* node, ISceneNode* parent) {
    IN("+ReadNode");
    ISceneNode* c = parent;

    // Read transformations
    Matrix4   _m = node->getTransformationMatrix();
    if (!_m.isIdentiy()) {
        TransformationNode* tn = new TransformationNode();
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
    }

    // Read instance geometry
    for (unsigned int i = 0; i < node->getInstanceGeometries().getCount(); i++) {
        ReadInstanceGeometry(node->getInstanceGeometries()[i], c);
    }

    // Read sub-nodes
    for (unsigned int i = 0; i < node->getInstanceNodes().getCount(); i++) {
        ReadNode(LookupNode(node->getInstanceNodes()[i]->getInstanciatedObjectId()), c);
    }
    for (unsigned int i = 0; i < node->getChildNodes().getCount(); i++) {
        ReadNode(node->getChildNodes()[i], c);
    }
    OUT();
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
    //logger.info << "start loading" << logger.end;
}

/** Remove all objects that don't have an object. Deletes unused visual scenes.*/
void ColladaResource::finish() {
    //logger.info << "finish loading" << logger.end;
}

/** When this method is called, the writer must write the global document asset.
	@return The writer should return true, if writing succeeded, false otherwise.*/
bool ColladaResource::writeGlobalAsset ( const COLLADAFW::FileInfo* asset ) {
    IN("writeGlobalAsset");
    coord[0] = 0;
    coord[1] = 1;
    coord[2] = 2;
    switch (asset->getUpAxisType()) {
    case FileInfo::X_UP:
        coord[0] = 1;
        coord[1] = 0;
        coord[2] = 2;
        break;
    case FileInfo::Z_UP:
        coord[0] = 0;
        coord[1] = 2;
        coord[2] = 1;
        break;
    case FileInfo::NONE:
    default:
        Warning("No up-axis defined. Assuming y is up");
    };
    OUT();
    return true;
}

/** Writes the entire visual scene.
    @return True on succeeded, false otherwise.*/
bool ColladaResource::writeVisualScene ( const COLLADAFW::VisualScene* visualScene ) {
    IN("Visual Scene");
    this->visualScene = new COLLADAFW::VisualScene(*visualScene);
    OUT();
    return true;
}

/** Writes the scene.
	@return True on succeeded, false otherwise.*/
bool ColladaResource::writeScene ( const COLLADAFW::Scene* scene ) {
    IN("writeScene");
    OUT();
    return true;
}

/** Handles all nodes in the library nodes.
	@return True on succeeded, false otherwise.*/
bool ColladaResource::writeLibraryNodes( const COLLADAFW::LibraryNodes* libraryNodes ) {
    IN("writeLibraryNodes");
    for (unsigned int i = 0; i < libraryNodes->getNodes().getCount(); i++) {
        Node* n = libraryNodes->getNodes()[i];
        nodes[n->getUniqueId()] = n;
    }
    OUT();
    return true;
}

/** Writes the geometry.
	@return True on succeeded, false otherwise.*/
bool ColladaResource::writeGeometry ( const COLLADAFW::Geometry* geometry ) {
    IN("Geometry");
    geometries[geometry->getUniqueId()] = ReadGeometry(geometry);
    OUT();
    return true;
}

/** Writes the material.
	@return True on succeeded, false otherwise.*/
bool ColladaResource::writeMaterial( const COLLADAFW::Material* material ) {
    IN("writeMaterial");
    materials[material->getUniqueId()] = new COLLADAFW::Material(*material);
    OUT();
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
    IN("writeEffect");

    CommonEffectPointerArray cepa = effect->getCommonEffects();
    if (cepa.getCount() == 0) {
        Warning("Effect is not a common effect.");
        return true;
    }
    EffectCommon* ce = cepa[0];
    MaterialPtr m = MaterialPtr(new Material());

    ExtractColor(ce->getAmbient(), m->ambient);
    ExtractColor(ce->getDiffuse(), m->diffuse);
    ExtractColor(ce->getSpecular(), m->specular);
    ExtractColor(ce->getEmission(), m->emission);
    ExtractFloatAttribute(ce->getShininess(), m->shininess);

    // logger.info << "ambient: " << m->ambient << logger.end;
    // logger.info << "diffuse: " << m->diffuse << logger.end;
    // logger.info << "specular: " << m->specular << logger.end;
    
    if (ce->getSamplerPointerArray().getCount() > 0) {
        // logger.info << "sampler is present" << logger.end;
        Sampler* s = ce->getSamplerPointerArray()[0];
        if (s->getSamplerType() != Sampler::SAMPLER_TYPE_2D) {
            Warning("Unsupported texture sampling type");
        }
        else {
            m->AddTexture(LookupImage(s->getSourceImage()));
        }
    }
    
    effects[effect->getUniqueId()] = m;
    OUT();
    // @todo: handle shader type
    return true;
}

    /** Writes the camera.
		@return True on succeeded, false otherwise.*/
bool ColladaResource::writeCamera( const COLLADAFW::Camera* camera ) {
    IN("writeCamera");
    OUT();
    return true;
}

    /** Writes the image.
		@return True on succeeded, false otherwise.*/
bool ColladaResource::writeImage( const COLLADAFW::Image* image ) {
    IN("writeImage");
    if (image->getSourceType() != Image::SOURCE_TYPE_URI) 
        Error("Unsupported image source type.");
    const URI& uri = image->getImageURI();
    try {
        ITextureResourcePtr texr = ResourceManager<ITextureResource>::Create(resource_dir + uri.getURIString());
        images[image->getUniqueId()] = texr;

    } catch (ResourceException e) {
        Warning("Error loading texture with message: " + string(e.what()));
    }
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
