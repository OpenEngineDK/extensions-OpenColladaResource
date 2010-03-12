// OpenCollada Model resource.
// -------------------------------------------------------------------
// Copyright (C) 2007 OpenEngine.dk (See AUTHORS) 
// Modified by Anders Bach Nielsen <abachn@daimi.au.dk> - 21. Nov 2007
// 
// This program is free software; It is covered by the GNU General 
// Public License version 2 or any later version. 
// See the GNU General Public License for more details (see LICENSE). 
//--------------------------------------------------------------------

#ifndef _OPEN_COLLADA_RESOURCE_H_
#define _OPEN_COLLADA_RESOURCE_H_

#include <Resources/IModelResource.h>
#include <Resources/IResourcePlugin.h>

#include <Math/Vector.h>

#include <string>
#include <map>
#include <list>

// OpenCollada stuff
#include <COLLADAFWIWriter.h>
#include <COLLADAFWUniqueId.h>
#include <COLLADAFWMeshVertexData.h>
#include <COLLADAFWInstanceGeometry.h>
#include <COLLADAFWMaterialBinding.h>

#include <Geometry/Mesh.h>

//forward declarations
namespace COLLADAFW {
    class VisualScene;
    class ColorOrTexture;
    class FloatOrParam;
    class Node;
    class Transformation;
    class Geometry;
}


namespace OpenEngine {
    //forward declarations
    namespace Geometry {
        class Material;
        class Mesh;
    }
    namespace Scene {
        class ISceneNode;
    }
namespace Resources {
    class ITexture2D;
    typedef boost::shared_ptr<ITexture2D> ITexture2DPtr;

using std::string;
using std::map;
using std::list;

using Geometry::Material;
using Geometry::Mesh;
using Math::Vector;
using Scene::ISceneNode;

typedef boost::shared_ptr<Material> MaterialPtr;

/**
 * OpenCollada resource.
 *
 * @class ColladaResource ColladaResource.h "ColladaResource.h"
 */
class ColladaResource : public IModelResource, public COLLADAFW::IWriter {
private: 
    class GeoPrimitive {
    public:
        COLLADAFW::MaterialId mId;
        Mesh* prim;
        GeoPrimitive(COLLADAFW::MaterialId mId, Mesh* prim) : mId(mId), prim(prim) {};
        virtual ~GeoPrimitive() { delete prim; }
    };
    typedef list<GeoPrimitive*> GeoPrimitives;
    string file, resource_dir, space;
    ISceneNode* root;
    COLLADAFW::VisualScene* visualScene;
    unsigned int coord[3]; // up-index correction vector
    
    // resource maps
    map<COLLADAFW::UniqueId, MaterialPtr>                effects;
    map<COLLADAFW::UniqueId, ITexture2DPtr>              images;
    map<COLLADAFW::UniqueId, GeoPrimitives*>             geometries;
    map<COLLADAFW::UniqueId, const COLLADAFW::Material*> materials;
    map<COLLADAFW::UniqueId, COLLADAFW::Node*>           nodes;

    // utility methods
    inline void Error(string msg);
    inline void Warning(string msg);

    inline COLLADAFW::Node* LookupNode(COLLADAFW::UniqueId id);
    inline GeoPrimitives*   LookupGeometry(COLLADAFW::UniqueId id);
    inline MaterialPtr      LookupMaterial(COLLADAFW::UniqueId id);
    inline ITexture2DPtr    LookupImage(COLLADAFW::UniqueId id);

    inline map<COLLADAFW::MaterialId, COLLADAFW::UniqueId> ExtractMaterialBindingMap(COLLADAFW::MaterialBindingArray& mbs);
    inline bool ExtractColor(COLLADAFW::ColorOrTexture& cot, Vector<4,float>& dest);
    inline bool ExtractFloatAttribute(COLLADAFW::FloatOrParam& fop, float& dest);
    inline bool ExtractFloatArray(COLLADAFW::MeshVertexData& d, float** dest);
    inline COLLADAFW::MeshVertexData::InputInfos* ExtractInputInfos(COLLADAFW::MeshVertexData& d);

    inline void           ReadNode(COLLADAFW::Node* node, ISceneNode* parent);
    inline void           ReadInstanceGeometry(COLLADAFW::InstanceGeometry* g, ISceneNode* parent);
    inline GeoPrimitives* ReadGeometry(const COLLADAFW::Geometry* g);
    
public:
    ColladaResource(string file);
    virtual ~ColladaResource();
    void Load();
    void Unload();
    ISceneNode* GetSceneNode();

    // IWriter CALLBACK Methods
    void cancel(const string& errorMessage);
    void start();
    void finish();

    bool writeGlobalAsset ( const COLLADAFW::FileInfo* asset );
    bool writeVisualScene ( const COLLADAFW::VisualScene* visualScene );
    bool writeScene ( const COLLADAFW::Scene* scene );
    bool writeLibraryNodes( const COLLADAFW::LibraryNodes* libraryNodes );
    bool writeGeometry ( const COLLADAFW::Geometry* geometry );
    bool writeMaterial( const COLLADAFW::Material* material );
    bool writeEffect( const COLLADAFW::Effect* effect );
    bool writeCamera( const COLLADAFW::Camera* camera );
    bool writeImage( const COLLADAFW::Image* image );
    bool writeLight( const COLLADAFW::Light* light );
    bool writeAnimation( const COLLADAFW::Animation* animation );
    bool writeAnimationList( const COLLADAFW::AnimationList* animationList );
    bool writeSkinControllerData( const COLLADAFW::SkinControllerData* skinControllerData );
    bool writeController( const COLLADAFW::Controller* Controller );
    bool writeFormulas( const COLLADAFW::Formulas* formulas );
    bool writeKinematicsScene( const COLLADAFW::KinematicsScene* kinematicsScene );
};

/**
 * Collada resource plug-in.
 *
 * @class ColladaPlugin ColladaResource.h "ColladaResource.h"
 */
class ColladaPlugin : public IResourcePlugin<IModelResource> {
public:
	ColladaPlugin();
    IModelResourcePtr CreateResource(string file);
};

} // NS Resources
} // NS OpenEngine

#endif // _OPEN_COLLADA_RESOURCE_H_
