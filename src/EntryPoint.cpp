#include <cassert>
#include <fbxsdk/fileio/fbxiosettings.h>
#include <fstream>
#include <iostream>
#include <fbxsdk.h>
#include <vector>

const int MDXM_IDENT = ('M' << 24) | ('G' << 16) | ('L' << 8) | '2';
const int MDXM_VERSION = 6;
const int MAX_QPATH = 64;

struct mdxmHeader_t
{
	int ident;
	int version;

	char name[MAX_QPATH];
	char animName[MAX_QPATH];

	int animIndex;
	int numBones;

	int numLODs;
	int ofsLODs;

	int numSurfaces;
	int ofsSurfHierarchy;

	int ofsEnd;
};

struct mdxmSurfHierarchy_t
{
	char name[MAX_QPATH];
	unsigned int flags;
	char shader[MAX_QPATH];
	int shaderIndex;
	int parentIndex;
	int numChildren;
	int childIndex[1];
};

struct mdxmSurface_t
{
	int ident;
	int thisSurfaceIndex;
	int ofsHeader;
	
	int numVerts;
	int ofsVerts;

	int numTriangles;
	int ofsTriangles;

	int numBoneReferences;
	int ofsBoneReferences;

	int ofsEnd;
};

struct mdxmTriangle_t
{
	int indexes[3];
};

struct mdxmVertex_t
{
	float normal[3];
	float position[3];
	unsigned int numWeightsAndBoneIndexes;
	unsigned char boneWeightings[4];
};

struct mdxmVertexTexcoord_t
{
	float st[2];
};

struct GLMSurfaceHierarchy
{
	const FbxNode *node;
	mdxmSurfHierarchy_t *metadata;
	std::size_t metadataSize;
};

struct GLMSurface
{
	mdxmSurface_t metadata;
	std::vector<mdxmTriangle_t> triangles;
	std::vector<mdxmVertex_t> vertices;
	std::vector<mdxmVertexTexcoord_t> texcoords;
};

struct ModelDetailData
{
	int lod;
	std::vector<GLMSurface> surfaces;
};

void PrintNodeAttribute ( const FbxNodeAttribute *attribute )
{
	switch ( attribute->GetAttributeType() )
	{
		case FbxNodeAttribute::eMesh:
		{
			std::cout << "There's a mesh!\n";

			const FbxMesh *mesh = static_cast<const FbxMesh *>(attribute);
			std::cout << "It has " << mesh->GetControlPointsCount() << " vertices.\n";
			std::cout << "It has " << mesh->mPolygons.Size() << " polygons\n";
			std::cout << "It has " << mesh->GetPolygonVertexCount() << " indices\n";

			if ( mesh->IsTriangleMesh() ) std::cout << "It has only triangles! :D\n";
			break;
		}

		default:
			std::cout << "There's no mesh!\n";
	}
}

void PrintNodeAttributes ( const FbxNode *node )
{
	// Default attribute first
	const FbxNodeAttribute *defaultAttribute = node->GetNodeAttribute();
	if ( defaultAttribute == NULL )
	{
		std::cout << "No default attribute for " << node->GetName() << '\n';
		return;
	}

	PrintNodeAttribute (defaultAttribute);

	for ( int i = 0; i < node->GetNodeAttributeCount(); i++ )
	{
		PrintNodeAttribute (node->GetNodeAttributeByIndex (i));
	}
}

const FbxMesh *GetFBXMesh ( const FbxNodeAttribute& attribute )
{
	switch ( attribute.GetAttributeType() )
	{
		case FbxNodeAttribute::eMesh:
			return static_cast<const FbxMesh *>(&attribute);
		default:
			return NULL;
	}

}

const FbxMesh *GetFBXMesh ( const FbxNode& node, int lod )
{
	// Default attribute first
	if ( lod == 0 )
	{
		const FbxNodeAttribute *defaultAttribute = node.GetNodeAttribute();
		if ( defaultAttribute != NULL )
		{
			return GetFBXMesh (*defaultAttribute);
		}

		return NULL;
	}

	const FbxMesh *mesh = NULL;
	for ( int i = 0; i < node.GetNodeAttributeCount(); i++ )
	{
		const FbxNodeAttribute *attribute = node.GetNodeAttributeByIndex (i);
		const FbxMesh *lodMesh = GetFBXMesh (*attribute);
		std::cout << "Trying attribute " << i << '\n';
		if ( lodMesh != NULL )
		{
			std::cout << "Found a mesh!\n";
			if ( lod == 0 )
			{
				break;
			}

			mesh = lodMesh;
			lod--;
		}
	}

	if ( lod > 0 )
	{
		return NULL;
	}

	return mesh;
}

void PrintSceneGraph ( const FbxNode *node, int level )
{
	for ( int i = 0; i < level; i++ )
	{
		std::cout << '.';
	}

	std::cout << node->GetName() << '\n';
	PrintNodeAttributes (node);

	for ( int i = 0; i < node->GetChildCount(); i++ )
	{
		PrintSceneGraph (node->GetChild (i), level + 1);
	}
}

const FbxNode *GetRootMesh ( const FbxNode& root )
{
	for ( int i = 0, count = root.GetChildCount(); i < count; i++ )
	{
		const FbxNode *child = root.GetChild (i);
		const FbxNodeAttribute *defaultAttribute = child->GetNodeAttribute();
		if ( defaultAttribute == NULL )
		{
			continue;
		}

		if ( defaultAttribute->GetAttributeType() == FbxNodeAttribute::eMesh )
		{
			return child;
		}
	}

	return NULL;
}

typedef std::vector<GLMSurfaceHierarchy> SurfaceHierarchyList;

int AddToHierarchy (
	int parentIndex,
	const FbxNode& node,
	SurfaceHierarchyList& hierarchy )
{
	int numChildren = node.GetChildCount();
	std::size_t hierarchySize = offsetof (mdxmSurfHierarchy_t, childIndex[numChildren]);
	char *data = new char[hierarchySize];
	mdxmSurfHierarchy_t *hierarchyNode = reinterpret_cast<mdxmSurfHierarchy_t *>(data);

	int thisIndex = hierarchy.size();

	if ( strncmp (node.GetName(), "bolt_", 5) == 0 )
	{
		hierarchyNode->name[0] = '*';
		strcpy (&hierarchyNode->name[1], node.GetName() + 5);
	}
	else
	{
		strcpy (hierarchyNode->name, node.GetName());
	}
	hierarchyNode->flags = 0;
	hierarchyNode->shader[0] = '\0';
	hierarchyNode->shaderIndex = 0;
	hierarchyNode->parentIndex = parentIndex;
	hierarchyNode->numChildren = numChildren;

	GLMSurfaceHierarchy surface;
	surface.metadata = hierarchyNode;
	surface.metadataSize = hierarchySize;
	surface.node = &node;

	hierarchy.push_back (surface);

	return thisIndex;
}

int CreateSurfaceHierarchy ( const FbxNode& node, SurfaceHierarchyList& hierarchyList, int parentIndex )
{
	int index = AddToHierarchy (parentIndex, node, hierarchyList);
	std::vector<int> childIndices;
	childIndices.reserve (node.GetChildCount());

	for ( int i = 0; i < node.GetChildCount(); i++ )
	{
		const FbxNode *child = node.GetChild (i);
		int childIndex = CreateSurfaceHierarchy (*child, hierarchyList, index);

		childIndices.push_back (childIndex);
	}

	for ( int i = 0; i < hierarchyList[index].metadata->numChildren; i++ )
	{
		hierarchyList[index].metadata->childIndex[i] = childIndices[i];
	}

	return index;
}

SurfaceHierarchyList CreateSurfaceHierarchy ( const FbxNode& root )
{
	SurfaceHierarchyList hierarchyList;

	CreateSurfaceHierarchy (root, hierarchyList, -1);

	return hierarchyList;
}

int GetLodCount ( const FbxNode& rootMesh )
{
	int lod = 0;
	while ( GetFBXMesh (rootMesh, lod) != NULL )
	{
		lod++;
	}

	return lod;
}

bool AllModelsHaveSameLoDs ( const SurfaceHierarchyList& hierarchy, int lod )
{
	for ( int i = 1; i < lod; i++ )
	{
		if ( GetLodCount (*hierarchy[i].node) != lod )
		{
			return false;
		}
	}

	return true;
}

std::size_t CalculateSurfaceSize ( const mdxmSurface_t& surface );
ModelDetailData CreateModelLod ( const SurfaceHierarchyList& hierarchy, int lod )
{
	ModelDetailData data;
	data.lod = lod;
	data.surfaces.resize (hierarchy.size());
	for ( int i = 0; i < hierarchy.size(); i++ )
	{
		const FbxNode& fbxNode = *hierarchy[i].node;
		const FbxMesh& mesh = *GetFBXMesh (fbxNode, lod);
		const mdxmSurfHierarchy_t& surfaceHierarchy = *hierarchy[i].metadata;
		mdxmSurface_t& metadata = data.surfaces[i].metadata;

		// Fill in header offset when writing.
		metadata.ident = 0;
		metadata.thisSurfaceIndex = i;
		metadata.ofsTriangles = sizeof (mdxmSurface_t);
		metadata.numTriangles = mesh.GetPolygonCount();
		metadata.numVerts = mesh.GetControlPointsCount();
		metadata.ofsVerts = metadata.ofsTriangles + sizeof (mdxmTriangle_t) * metadata.numTriangles;
		metadata.numBoneReferences = 0;
		metadata.ofsBoneReferences = 0;
		metadata.ofsEnd = CalculateSurfaceSize (metadata);
		
		// Triangles
		std::vector<mdxmTriangle_t>& triangles = data.surfaces[i].triangles;
		triangles.resize (metadata.numTriangles);

		const int *indices = mesh.GetPolygonVertices();
		for ( int j = 0, count = triangles.size(); j < count; j++, indices += 3 )
		{
			triangles[j].indexes[0] = indices[0];
			triangles[j].indexes[1] = indices[1];
			triangles[j].indexes[2] = indices[2];
		}

		// Vertices
		const FbxVector4 *positions = mesh.GetControlPoints();
		const FbxLayer *layer0 = mesh.GetLayer (0);
		if ( layer0 == NULL )
		{
			std::cerr << "No layer 0 for mesh " << i << " (" << surfaceHierarchy.name << ")\n";
		}

		const FbxLayerElementNormal *normals = layer0->GetNormals();
		if ( normals == NULL )
		{
			std::cerr << "No normals for mesh " << i << " (" << surfaceHierarchy.name << ")\n";
		}

		std::vector<mdxmVertex_t>& vertices = data.surfaces[i].vertices;
		vertices.resize (metadata.numVerts);
		for ( int i = 0, count = metadata.numVerts; i < count; i++ )
		{
			mdxmVertex_t& vertex = vertices[i];
			vertex.position[0] = static_cast<float>(positions[i][0]);
			vertex.position[1] = static_cast<float>(positions[i][1]);
			vertex.position[2] = static_cast<float>(positions[i][2]);

			FbxVector4 normal = normals->GetDirectArray().GetAt (i);
			vertex.normal[0] = normal[0];
			vertex.normal[1] = normal[1];
			vertex.normal[2] = normal[2];

			vertex.boneWeightings[0] = 0;
			vertex.boneWeightings[1] = 0;
			vertex.boneWeightings[2] = 0;
			vertex.boneWeightings[3] = 0;
			vertex.numWeightsAndBoneIndexes = 0;
		}

		// Texcoords
		const FbxLayerElementUV *uvs = layer0->GetUVs();

		std::vector<mdxmVertexTexcoord_t>& texcoords = data.surfaces[i].texcoords;
		texcoords.resize (metadata.numVerts);
		if ( uvs == NULL )
		{
			// Not all surfaces have/need texcoords. For example bolts, and tags.
			std::memset (&texcoords[0], 0, sizeof (mdxmVertexTexcoord_t));
		}
		else
		{
			for ( int i = 0, count = metadata.numVerts; i < count; i++ )
			{
				mdxmVertexTexcoord_t& texcoord = texcoords[i];

				FbxVector2 tc = uvs->GetDirectArray().GetAt (i);
				texcoord.st[0] = tc[0];
				texcoord.st[1] = tc[1];
			}
		}
	}

	return data;
}

std::vector<ModelDetailData> GetModelData ( const SurfaceHierarchyList& hierarchy )
{
	std::vector<ModelDetailData> detailData;
	if ( hierarchy.size() == 0 )
	{
		return detailData;
	}

	const FbxNode& root = *hierarchy[0].node;
	int lodCount = GetLodCount (root);

	if ( !AllModelsHaveSameLoDs (hierarchy, lodCount) )
	{
		std::cerr << "Some surfaces have more level of details than others.\n";
		return detailData;
	}

	detailData.reserve (lodCount);
	std::cout << "There are " << lodCount << " LoDs in root mesh\n";
	for ( int i = 0; i < lodCount; i++ )
	{
		detailData.push_back (CreateModelLod (hierarchy, i));
	}

	return detailData;
}

std::size_t CalculateHeaderSize()
{
	return sizeof (mdxmHeader_t);
}

std::size_t CalculateSurfaceHierarchySize ( const SurfaceHierarchyList& hierarchy )
{
	std::size_t filesize = 0;

	// Hierarchy offsets
	filesize += hierarchy.size() * sizeof (int);

	// Surface hierarchy data
	for ( int i = 0, count = hierarchy.size(); i < count; i++ )
	{
		filesize += hierarchy[i].metadataSize;
	}

	return filesize;
}

std::size_t CalculateSurfaceSize ( const mdxmSurface_t& surface )
{
	std::size_t filesize = 0;

	// Surface data
	filesize += sizeof (mdxmSurface_t);

	// Triangles
	filesize += sizeof (mdxmTriangle_t) * surface.numTriangles;

	// Vertices
	filesize += sizeof (mdxmVertex_t) * surface.numVerts;

	// Texcoords
	filesize += sizeof (mdxmVertexTexcoord_t) * surface.numVerts;

	return filesize;
}

std::size_t CalculateLODSize (
		const std::vector<ModelDetailData>& modelData,
		int lod )
{
	std::size_t filesize = 0;

	// LOD data
	const std::vector<GLMSurface>& surfaces = modelData[lod].surfaces;
	
	// Next LOD offset
	filesize += sizeof (int);

	// Surface offsets
	filesize += sizeof (int) * modelData[lod].surfaces.size();

	for ( int i = 0, count = modelData[lod].surfaces.size(); i < count; i++ )
	{
		const GLMSurface& surface = surfaces[i];

		filesize += CalculateSurfaceSize (surface.metadata);
	}

	return filesize;

}

std::size_t CalculateLODSize (
		const std::vector<ModelDetailData>& modelData )
{
	std::size_t filesize = 0;

	for ( int lod = 0, count = modelData.size(); lod < count; lod++ )
	{
		filesize += CalculateLODSize (modelData, lod);
	}

	return filesize;
}

std::size_t CalculateGLMFileSize (
		const SurfaceHierarchyList& hierarchy,
		const std::vector<ModelDetailData>& modelData )
{
	return CalculateHeaderSize () + 
			CalculateSurfaceHierarchySize (hierarchy) +
			CalculateLODSize (modelData);
}

void MakeGLMFile ( const FbxNode& root )
{
	const FbxNode *meshRoot = GetRootMesh (root);

	// Create surface hierarchy
	SurfaceHierarchyList surfaceHierarchy (CreateSurfaceHierarchy (*meshRoot));

	// Create surface data
	std::vector<ModelDetailData> modelDetails (GetModelData (surfaceHierarchy));

	std::size_t filesize = CalculateGLMFileSize (surfaceHierarchy, modelDetails);
	
	// Write data
	std::vector<char> buffer (filesize, '\0');

	int hierarchyOffsetsBase = sizeof (mdxmHeader_t);
	int hierarchyBase = hierarchyOffsetsBase + sizeof (int) * surfaceHierarchy.size();
	int lodBase = hierarchyOffsetsBase + CalculateSurfaceHierarchySize (surfaceHierarchy);

	mdxmHeader_t *header = reinterpret_cast<mdxmHeader_t *>(&buffer[0]);
	header->ident = MDXM_IDENT;
	header->version = MDXM_VERSION;
	strcpy (header->name, "test.glm");
	strcpy (header->animName, "models/players/_humanoid/_humanoid");
	header->animIndex = -1;
	header->numBones = 53;
	header->numLODs = static_cast<int>(modelDetails.size());
	header->ofsLODs = lodBase;
	header->numSurfaces = static_cast<int>(surfaceHierarchy.size());
	header->ofsSurfHierarchy = hierarchyBase;
	header->ofsEnd = filesize;

	// Surface offsets and hierarchy
	char *hierarchyOffsetsBasePtr = &buffer[hierarchyOffsetsBase];
	int *hierarchyOffsets = reinterpret_cast<int *>(hierarchyOffsetsBasePtr);

	// Base for the offset is at the start of this array of offsets.
	int hierarchyOffset = sizeof (int) * header->numSurfaces;
	mdxmSurfHierarchy_t *hierarchy = reinterpret_cast<mdxmSurfHierarchy_t *>(&hierarchyOffsetsBasePtr[hierarchyOffset]);

	for ( int i = 0; i < header->numSurfaces; i++ )
	{
		const GLMSurfaceHierarchy& surface = surfaceHierarchy[i];

		hierarchyOffsets[i] = hierarchyOffset;

		std::memcpy (hierarchy, surface.metadata, surface.metadataSize);

		hierarchyOffset += surface.metadataSize;
		hierarchy = reinterpret_cast<mdxmSurfHierarchy_t *>(&hierarchyOffsetsBasePtr[hierarchyOffset]);
	}

	// LODs
	char *lodBasePtr = reinterpret_cast<char *>(&buffer[lodBase]);
	int lodOffset = 0;
	for ( int lod = 0; lod < header->numLODs; lod++ )
	{
		const ModelDetailData& modelDetail = modelDetails[lod];

		// Write offset to next header
		int lodOffset = CalculateLODSize (modelDetails, lod);
		int *nextLOD = reinterpret_cast<int *>(lodBasePtr);
		*nextLOD = lodOffset;

		int *surfaceOffsets = reinterpret_cast<int *>(lodBasePtr + sizeof (int));
		char *surfaceData = reinterpret_cast<char *>(surfaceOffsets);
		int surfaceOffset = sizeof (int) * header->numSurfaces;

		for ( int i = 0; i < header->numSurfaces; i++ )
		{
			const GLMSurface& glmSurface = modelDetail.surfaces[i];

			surfaceOffsets[i] = surfaceOffset;

			mdxmSurface_t *surface = reinterpret_cast<mdxmSurface_t *>(&surfaceData[surfaceOffset]);
			*surface = glmSurface.metadata;
			surface->ofsHeader = &buffer[0] - reinterpret_cast<char *>(surface);
			surfaceOffset += sizeof (mdxmSurface_t);

			mdxmTriangle_t *triangles = reinterpret_cast<mdxmTriangle_t *>(&surfaceData[surfaceOffset]);
			std::memcpy (triangles, glmSurface.triangles.data(), sizeof (mdxmTriangle_t) * surface->numTriangles);
			surfaceOffset += sizeof (mdxmTriangle_t) * surface->numTriangles;

			mdxmVertex_t *vertices = reinterpret_cast<mdxmVertex_t *>(&surfaceData[surfaceOffset]);
			std::memcpy (vertices, glmSurface.vertices.data(), sizeof (mdxmVertex_t) * surface->numVerts);
			surfaceOffset += sizeof (mdxmVertex_t) * surface->numVerts;

			mdxmVertexTexcoord_t *texcoords = reinterpret_cast<mdxmVertexTexcoord_t *>(&surfaceData[surfaceOffset]);
			std::memcpy (texcoords, glmSurface.texcoords.data(), sizeof (mdxmVertexTexcoord_t) * surface->numVerts);
			surfaceOffset += sizeof (mdxmVertexTexcoord_t) * surface->numVerts;
		}

		lodBasePtr += lodOffset;
	}

	std::ofstream file ("test.glm", std::ios::binary);
	if ( !file )
	{
		std::cerr << "Failed to create file test.glm\n";
	}

	file.write (buffer.data(), buffer.size());
}

int main()
{
	std::cout << "Hello World\n";
	FbxManager *fbxManager = FbxManager::Create();

	FbxIOSettings *ios = FbxIOSettings::Create (fbxManager, IOSROOT);
	fbxManager->SetIOSettings (ios);

	FbxImporter *importer = FbxImporter::Create (fbxManager, "");
	const char *modelPath = "model.fbx";

	bool importStatus = importer->Initialize (modelPath, -1, fbxManager->GetIOSettings());
	if ( !importStatus )
	{
		std::cout << "Call to FbxImporter::Initialize() failed.\n";
		//std::cout << "Error returned: " << importer->GetLastErrorString() << '\n';

		return 1;
	}

	FbxScene *scene = FbxScene::Create (fbxManager, "model");
	importer->Import (scene);

	if ( scene->GetCharacterCount() > 0 )
	{
		std::cout << "There are " << scene->GetCharacterCount() << " characters\n";
	}
	else
	{
		std::cout << "No character!\n";
	}

	FbxNode *root = scene->GetRootNode();
	if ( root == NULL )
	{
		std::cout << "No root node :(\n";
		return 1;
	}

	//PrintSceneGraph (root, 0);
	MakeGLMFile (*root);

	importer->Destroy();
}
