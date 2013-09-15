#include <cassert>
#include <fbxsdk/fileio/fbxiosettings.h>
#include <fstream>
#include <iostream>
#include <fbxsdk.h>
#include <map>
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

const int G2SURFACE_BOLT		= 0x1;
const int G2SURFACE_HIDDEN	= 0x2;

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
	FbxNode *node;
	mdxmSurfHierarchy_t *metadata;
	std::size_t metadataSize;
};

struct GLMSurface
{
	mdxmSurface_t metadata;
	std::vector<mdxmTriangle_t> triangles;
	std::vector<mdxmVertex_t> vertices;
	std::vector<mdxmVertexTexcoord_t> texcoords;
	std::vector<int> boneReferences;
};

struct ModelDetailData
{
	int lod;
	std::vector<GLMSurface> surfaces;
};

void PrintNodeAttribute ( FbxNodeAttribute *attribute )
{
	switch ( attribute->GetAttributeType() )
	{
		case FbxNodeAttribute::eMesh:
		{
			std::cout << "There's a mesh!\n";

			FbxMesh *mesh = static_cast<FbxMesh *>(attribute);
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

void PrintNodeAttributes ( FbxNode *node )
{
	// Default attribute first
	FbxNodeAttribute *defaultAttribute = node->GetNodeAttribute();
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

FbxMesh *GetFBXMesh ( FbxNodeAttribute& attribute )
{
	switch ( attribute.GetAttributeType() )
	{
		case FbxNodeAttribute::eMesh:
			return static_cast<FbxMesh *>(&attribute);
		default:
			return NULL;
	}

}

FbxMesh *GetFBXMesh ( FbxNode& node, int lod )
{
	// Default attribute first
	if ( lod == 0 )
	{
		FbxNodeAttribute *defaultAttribute = node.GetNodeAttribute();
		if ( defaultAttribute != NULL )
		{
			return GetFBXMesh (*defaultAttribute);
		}

		return NULL;
	}

	FbxMesh *mesh = NULL;
	for ( int i = 0; i < node.GetNodeAttributeCount(); i++ )
	{
		FbxNodeAttribute *attribute = node.GetNodeAttributeByIndex (i);
		FbxMesh *lodMesh = GetFBXMesh (*attribute);
		if ( lodMesh != NULL )
		{
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

void PrintSceneGraph ( FbxNode *node, int level )
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

void PrintAllPropertyNames ( const FbxObject& object )
{
	FbxProperty prop = object.GetFirstProperty();
	do
	{
		std::cout << prop.GetName() << '\n';
		prop = object.GetNextProperty (prop);
	} while ( prop.IsValid() );

	std::cout << "\n\n";
}

FbxNode *GetRootMesh ( FbxNode& root )
{
	for ( int i = 0, count = root.GetChildCount(); i < count; i++ )
	{
		FbxNode *child = root.GetChild (i);
		FbxNodeAttribute *defaultAttribute = child->GetNodeAttribute();
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
	FbxNode& node,
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

	if ( hierarchyNode->name[0] == '*' )
	{
		hierarchyNode->flags |= G2SURFACE_BOLT;
	}

	if ( strstr (hierarchyNode->name, "_off") )
	{
		hierarchyNode->flags |= G2SURFACE_HIDDEN;
	}

	GLMSurfaceHierarchy surface;
	surface.metadata = hierarchyNode;
	surface.metadataSize = hierarchySize;
	surface.node = &node;

	hierarchy.push_back (surface);

	return thisIndex;
}

int CreateSurfaceHierarchy ( FbxNode& node, SurfaceHierarchyList& hierarchyList, int parentIndex )
{
	int index = AddToHierarchy (parentIndex, node, hierarchyList);
	std::vector<int> childIndices;
	childIndices.reserve (node.GetChildCount());

	for ( int i = 0; i < node.GetChildCount(); i++ )
	{
		FbxNode *child = node.GetChild (i);
		int childIndex = CreateSurfaceHierarchy (*child, hierarchyList, index);

		childIndices.push_back (childIndex);
	}

	for ( int i = 0; i < hierarchyList[index].metadata->numChildren; i++ )
	{
		hierarchyList[index].metadata->childIndex[i] = childIndices[i];
	}

	return index;
}

SurfaceHierarchyList CreateSurfaceHierarchy ( FbxNode& root )
{
	SurfaceHierarchyList hierarchyList;

	CreateSurfaceHierarchy (root, hierarchyList, -1);

	return hierarchyList;
}

int GetLodCount ( FbxNode& rootMesh )
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

struct VertexId
{
	VertexId() {}
	VertexId ( int positionId, int texcoordId )
		: positionId (positionId),
			texcoordId (texcoordId)
	{
	}

	int positionId;
	int texcoordId;

	bool operator < ( const VertexId& id ) const 
	{
		if ( positionId < id.positionId )
		{
			return true;
		}
		else if ( positionId > id.positionId )
		{
			return false;
		}
		else
		{
			return texcoordId < id.texcoordId;
		}
	}
};

struct Vertex
{
	Vertex(): filled (false) {}
	bool filled;
	mdxmVertex_t positionAndNormal;
	mdxmVertexTexcoord_t texcoord;
};

std::size_t CalculateSurfaceSize ( const mdxmSurface_t& surface );
ModelDetailData CreateModelLod ( FbxScene& scene, const SurfaceHierarchyList& hierarchy, int lod )
{
	ModelDetailData data;
	data.lod = lod;
	data.surfaces.resize (hierarchy.size());

	std::vector<int> newIndices;
	std::map<VertexId, int> uniqueVerticesMap;
	std::vector<Vertex> uniqueVertices;
	for ( int i = 0; i < hierarchy.size(); i++ )
	{
		FbxNode& fbxNode = *hierarchy[i].node;
		FbxMesh& mesh = *GetFBXMesh (fbxNode, lod);
		const mdxmSurfHierarchy_t& surfaceHierarchy = *hierarchy[i].metadata;
		mdxmSurface_t& metadata = data.surfaces[i].metadata;

		// Get the relevant vertex data from the mesh first.
		FbxAMatrix globalMatrix;
		globalMatrix = scene.GetEvaluator()->GetNodeGlobalTransform (&fbxNode);

		const int *indices = mesh.GetPolygonVertices();
		const FbxVector4 *positions = mesh.GetControlPoints();
		const FbxLayer *layer0 = mesh.GetLayer (0);
		if ( layer0 == NULL )
		{
			std::cerr << "No layer 0 for mesh " << i << " (" << surfaceHierarchy.name << ")\n";
		}

		const FbxLayerElementNormal *normalsLayer = layer0->GetNormals();
		if ( normalsLayer == NULL )
		{
			std::cerr << "No normals for mesh " << i << " (" << surfaceHierarchy.name << ")\n";
		}

		FbxLayerElementArrayTemplate<FbxVector4>& normals = normalsLayer->GetDirectArray();

		const FbxLayerElementUV *uvLayer = layer0->GetUVs();
		FbxLayerElementArrayTemplate<FbxVector2> *uvs = NULL;

		int numUniquePoints = mesh.GetControlPointsCount();
		int numPositions = mesh.GetControlPointsCount();
		if ( uvLayer != NULL )
		{
			if ( !mesh.GetTextureUV (&uvs) )
			{
				std::cerr << "Failed to retrieve texture UVs when UV layer exists.\n";
			}
			else
			{
				numUniquePoints = uvs->GetCount();
			}
		}

		int numTriangles = mesh.GetPolygonCount();
		int numVerts = numUniquePoints;

		std::vector<mdxmTriangle_t>& triangles = data.surfaces[i].triangles;
		triangles.resize (numTriangles);

		uniqueVertices.clear();
		if ( uniqueVertices.size() < numUniquePoints )
		{
			uniqueVertices.resize (numUniquePoints);
		}

		if ( uvs == NULL )
		{
			// The only surfaces with no UVs are the tags and the "stupidtriangle" surface.
			
			// Tags are oriented by taking the last vertex as the origin, and calculating
			// the different edge lengths. The mid-length edge is taken to be the Y axis,
			// and shortest to be X axis.
			int order[] = {1, 2, 0};
			for ( int j = 0, count = 3; j < count; j++ )
			{
				Vertex& vertex = uniqueVertices[j];
				if ( vertex.filled )
				{
					continue;
				}

				vertex.filled = true;

				FbxVector4 position = globalMatrix.MultT (positions[order[j]]);
				vertex.positionAndNormal.position[0] = static_cast<float>(position[0]);
				vertex.positionAndNormal.position[1] = static_cast<float>(position[1]);
				vertex.positionAndNormal.position[2] = static_cast<float>(position[2]);

				FbxVector4 normal = globalMatrix.MultT (normals[order[j]]);
				vertex.positionAndNormal.normal[0] = normal[0];
				vertex.positionAndNormal.normal[1] = normal[1];
				vertex.positionAndNormal.normal[2] = normal[2];

				vertex.positionAndNormal.numWeightsAndBoneIndexes = 0;
				vertex.positionAndNormal.boneWeightings[0] = 255;
				vertex.positionAndNormal.boneWeightings[1] = 0;
				vertex.positionAndNormal.boneWeightings[2] = 0;
				vertex.positionAndNormal.boneWeightings[3] = 0;

				// Not strictly necessary but makes things cleaner in the output file.
				vertex.texcoord.st[0] = 0.0f;
				vertex.texcoord.st[1] = 0.0f;
			}

			// Indices
			triangles[0].indexes[0] = 0;
			triangles[0].indexes[1] = 1;
			triangles[0].indexes[2] = 2;
		}
		else
		{
			assert (uvLayer->GetMappingMode() == FbxLayerElement::eByPolygonVertex &&
					uvLayer->GetReferenceMode() == FbxLayerElement::eIndexToDirect);

			FbxLayerElementArrayTemplate<int>& uvIndices = uvLayer->GetIndexArray();
			if ( numUniquePoints == numPositions )
			{
				for ( int j = 0, count = uvIndices.GetCount(); j < count; j++ )
				{
					Vertex& vertex = uniqueVertices[indices[j]];
					if ( vertex.filled )
					{
						continue;
					}

					vertex.filled = true;

					FbxVector4 position = globalMatrix.MultT (positions[indices[j]]);
					vertex.positionAndNormal.position[0] = static_cast<float>(position[0]);
					vertex.positionAndNormal.position[1] = static_cast<float>(position[1]);
					vertex.positionAndNormal.position[2] = static_cast<float>(position[2]);

					FbxVector4 normal = globalMatrix.MultT (normals[indices[j]]);
					vertex.positionAndNormal.normal[0] = normal[0];
					vertex.positionAndNormal.normal[1] = normal[1];
					vertex.positionAndNormal.normal[2] = normal[2];

					vertex.positionAndNormal.numWeightsAndBoneIndexes = 0;
					vertex.positionAndNormal.boneWeightings[0] = 255;
					vertex.positionAndNormal.boneWeightings[1] = 0;
					vertex.positionAndNormal.boneWeightings[2] = 0;
					vertex.positionAndNormal.boneWeightings[3] = 0;

					FbxVector2 tc = uvs->GetAt (uvIndices[j]);
					vertex.texcoord.st[0] = tc[0];
					vertex.texcoord.st[1] = 1.0f - tc[1];
				}

				// Triangles
				const int *index = indices;
				for ( int j = 0, count = triangles.size(); j < count; j++, index += 3 )
				{
					triangles[j].indexes[0] = index[0];
					triangles[j].indexes[1] = index[1];
					triangles[j].indexes[2] = index[2];
				}
			}
			else
			{
				uniqueVerticesMap.clear();
				int index = 0;
				for ( int tri = 0, k = 0; tri < numTriangles; tri++ )
				{
					int triangle[3];
					for ( int j = 0; j < 3; k++, j++ )
					{
						VertexId id (indices[k], uvIndices[k]);
						std::map<VertexId, int>::const_iterator it = uniqueVerticesMap.find (id);

						if ( it == uniqueVerticesMap.end() )
						{
							uniqueVerticesMap.insert (std::make_pair (id, index));
							Vertex& v = uniqueVertices[index];

							FbxVector4 position = globalMatrix.MultT (positions[id.positionId]);
							v.positionAndNormal.position[0] = static_cast<float>(position[0]);
							v.positionAndNormal.position[1] = static_cast<float>(position[1]);
							v.positionAndNormal.position[2] = static_cast<float>(position[2]);

							FbxVector4 normal = globalMatrix.MultT (normals[id.positionId]);
							v.positionAndNormal.normal[0] = normal[0];
							v.positionAndNormal.normal[1] = normal[1];
							v.positionAndNormal.normal[2] = normal[2];

							v.positionAndNormal.numWeightsAndBoneIndexes = 0;
							v.positionAndNormal.boneWeightings[0] = 255;
							v.positionAndNormal.boneWeightings[1] = 0;
							v.positionAndNormal.boneWeightings[2] = 0;
							v.positionAndNormal.boneWeightings[3] = 0;

							FbxVector2 tc = uvs->GetAt (id.texcoordId);
							v.texcoord.st[0] = tc[0];
							v.texcoord.st[1] = 1.0f - tc[1];

							triangle[j] = index;
							index++;
						}
						else
						{
							triangle[j] = it->second;
						}
					}

					triangles[tri].indexes[0] = triangle[0];
					triangles[tri].indexes[1] = triangle[1];
					triangles[tri].indexes[2] = triangle[2];
				}

				numVerts = index;
			}
		}

		std::vector<mdxmVertexTexcoord_t>& texcoords = data.surfaces[i].texcoords;
		texcoords.resize (numVerts);

		std::vector<mdxmVertex_t>& vertices = data.surfaces[i].vertices;
		vertices.resize (numVerts);
		for ( int j = 0; j < numVerts; j++ )
		{
			vertices[j] = uniqueVertices[j].positionAndNormal;
			texcoords[j] = uniqueVertices[j].texcoord;
		}

		// Bone references
		std::vector<int>& boneReferences = data.surfaces[i].boneReferences;
		boneReferences.resize (1);
		boneReferences[0] = 0;

		// Fill in header offset when writing.
		metadata.ident = 0;
		metadata.thisSurfaceIndex = i;
		metadata.ofsTriangles = sizeof (mdxmSurface_t);
		metadata.numTriangles = numTriangles;
		metadata.numVerts = numVerts;
		metadata.ofsVerts = metadata.ofsTriangles + sizeof (mdxmTriangle_t) * metadata.numTriangles;
		metadata.numBoneReferences = 1;

		int ofsBoneReferences = 0;
		ofsBoneReferences += sizeof (mdxmSurface_t);
		ofsBoneReferences += sizeof (mdxmTriangle_t) * metadata.numTriangles;
		ofsBoneReferences += sizeof (mdxmVertex_t) * metadata.numVerts;
		ofsBoneReferences += sizeof (mdxmVertexTexcoord_t) * metadata.numVerts;

		metadata.ofsBoneReferences = ofsBoneReferences;
		metadata.ofsEnd = CalculateSurfaceSize (metadata);
	}

	return data;
}

std::vector<ModelDetailData> GetModelData ( FbxScene& scene, const SurfaceHierarchyList& hierarchy )
{
	std::vector<ModelDetailData> detailData;
	if ( hierarchy.size() == 0 )
	{
		return detailData;
	}

	FbxNode& root = *hierarchy[0].node;
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
		detailData.push_back (CreateModelLod (scene, hierarchy, i));
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

	// Bone references!
	filesize += sizeof (int) * surface.numBoneReferences;

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

void MakeGLMFile ( FbxScene& scene, FbxNode& root )
{
	FbxNode *meshRoot = GetRootMesh (root);

	// Create surface hierarchy
	SurfaceHierarchyList surfaceHierarchy (CreateSurfaceHierarchy (*meshRoot));

	// Create surface data
	std::vector<ModelDetailData> modelDetails (GetModelData (scene, surfaceHierarchy));

	std::size_t filesize = CalculateGLMFileSize (surfaceHierarchy, modelDetails);
	
	// Write data
	std::vector<char> buffer (filesize, '\0');

	int hierarchyOffsetsBase = sizeof (mdxmHeader_t);
	int hierarchyBase = hierarchyOffsetsBase + sizeof (int) * surfaceHierarchy.size();
	int lodBase = hierarchyOffsetsBase + CalculateSurfaceHierarchySize (surfaceHierarchy);

	mdxmHeader_t *header = reinterpret_cast<mdxmHeader_t *>(&buffer[0]);
	header->ident = MDXM_IDENT;
	header->version = MDXM_VERSION;
	strcpy (header->name, "model.glm");
	strcpy (header->animName, "*default");
	header->animIndex = -1;
	header->numBones = 1;
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

			int *boneReferences = reinterpret_cast<int *>(&surfaceData[surfaceOffset]);
			std::memcpy (boneReferences, glmSurface.boneReferences.data(), sizeof (int) * surface->numBoneReferences);
			surfaceOffset += sizeof (int) * surface->numBoneReferences;
		}

		lodBasePtr += lodOffset;
	}

	for ( int i = 0; i < surfaceHierarchy.size(); i++ )
	{
		char *data = reinterpret_cast<char *>(surfaceHierarchy[i].metadata);
		delete [] data;
	}

	std::ofstream file ("model.glm", std::ios::binary);
	if ( !file )
	{
		std::cerr << "Failed to create file model.glm\n";
	}

	file.write (buffer.data(), buffer.size());
}

int main ( int argc, char *argv[] )
{
	std::cout << "Hello World\n";
	FbxManager *fbxManager = FbxManager::Create();

	FbxIOSettings *ios = FbxIOSettings::Create (fbxManager, IOSROOT);
	fbxManager->SetIOSettings (ios);

	FbxImporter *importer = FbxImporter::Create (fbxManager, "");
	const char *modelPath = "model.fbx";
	if ( argc > 1 )
	{
		modelPath = argv[1];
		std::cout << "Converting " << modelPath << '\n';
	}

	bool importStatus = importer->Initialize (modelPath, -1, fbxManager->GetIOSettings());
	if ( !importStatus )
	{
		std::cout << "Call to FbxImporter::Initialize() failed.\n";

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

	MakeGLMFile (*scene, *root);

	importer->Destroy();
}
