#include <algorithm>
#include <cassert>
#include <cctype>
#include <fbxsdk/fileio/fbxiosettings.h>
#include <fstream>
#include <iostream>
#include <fbxsdk.h>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "GLA.h"
#include "GLM.h"
#include "IO.h"

namespace
{

FbxMesh *GetFBXMesh ( FbxNodeAttribute& attribute )
{
	switch ( attribute.GetAttributeType() )
	{
	case FbxNodeAttribute::eMesh:
		return static_cast<FbxMesh *>(&attribute);
	default:
		return nullptr;
	}
}

FbxMesh *GetFBXMesh ( FbxNode& node )
{
	for ( int i = 0; i < node.GetNodeAttributeCount(); i++ )
	{
		FbxNodeAttribute *attribute = node.GetNodeAttributeByIndex (i);
		FbxMesh *lodMesh = GetFBXMesh (*attribute);

		if ( lodMesh != nullptr )
		{
			return lodMesh;
		}
	}

	return nullptr;
}

char *CopyString ( char *destination, const char *source, std::size_t destinationSize )
{
	std::strncpy (destination, source, destinationSize - 1);
	destination[destinationSize - 1] = '\0';

	return destination;
}

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
		CopyString (hierarchyNode->name + 1, node.GetName() + 5, sizeof (hierarchyNode->name) - 1);
	}
	else
	{
		CopyString (hierarchyNode->name, node.GetName(), sizeof (hierarchyNode->name));
	}
	std::transform (hierarchyNode->name, hierarchyNode->name + sizeof (hierarchyNode->name),
		hierarchyNode->name, []( char c ) { return std::tolower (c); });

	hierarchyNode->shader[0] = '\0';
	if (node.GetSrcObjectCount<FbxSurfacePhong>() > 0)
	{
		FbxSurfacePhong *material = static_cast<FbxSurfacePhong *>(node.GetSrcObject<FbxSurfacePhong>());
		if (material->Diffuse.GetSrcObjectCount<FbxFileTexture>() > 0)
		{
			CopyString (hierarchyNode->shader, material->GetInitialName(), sizeof (hierarchyNode->shader));
		}
	}

	hierarchyNode->flags = 0;
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

	hierarchy.push_back (surface);

	return thisIndex;
}

int CreateSurfaceHierarchy ( FbxNode& node, SurfaceHierarchyList& hierarchyList, int parentIndex )
{
	int index = 0;
	std::vector<int> childIndices;

	index = AddToHierarchy (parentIndex, node, hierarchyList);
	childIndices.reserve (node.GetChildCount());

	for ( int i = 0; i < node.GetChildCount(); i++ )
	{
		FbxNode *child = node.GetChild (i);
		int childIndex = CreateSurfaceHierarchy (*child, hierarchyList, index);

		childIndices.push_back (childIndex);
	}

	std::copy (childIndices.begin(), childIndices.end(),
		&hierarchyList[index].metadata->childIndex[0]);

	return index;
}

SurfaceHierarchyList CreateSurfaceHierarchy ( FbxNode& root )
{
	SurfaceHierarchyList hierarchyList;

	CreateSurfaceHierarchy (root, hierarchyList, -1);

	return hierarchyList;
}

struct VertexId
{
	VertexId() {}
	VertexId ( int positionId, int texcoordId )
		: positionId (positionId)
		, texcoordId (texcoordId)
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

void VectorNormalize ( float *v )
{
	float dot = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	if ( dot < 1e-6f )
	{
		return;
	}

	float len = 1.0f / std::sqrt (dot);
	v[0] *= len;
	v[1] *= len;
	v[2] *= len;
}

struct Weights
{
	std::string influencers[4];
	float weights[4];
	std::size_t count;
};

void CopyVertexData (
	Vertex& vertex,
	const FbxVector4& position,
	const FbxVector4& normal,
	const Weights& weights,
	const std::map<std::string, int>& boneIndexes )
{
	vertex.positionAndNormal.position[0] = static_cast<float>(position[0]);
	vertex.positionAndNormal.position[1] = static_cast<float>(position[1]);
	vertex.positionAndNormal.position[2] = static_cast<float>(position[2]);

	vertex.positionAndNormal.normal[0] = static_cast<float>(normal[0]);
	vertex.positionAndNormal.normal[1] = static_cast<float>(normal[1]);
	vertex.positionAndNormal.normal[2] = static_cast<float>(normal[2]);
	VectorNormalize (vertex.positionAndNormal.normal);

	vertex.positionAndNormal.numWeightsAndBoneIndexes = 0;

	if ( weights.count > 4 )
	{
		std::cerr << "Bad things happened.\n";
		return;
	}

	int bitOffset = 0;
	int overflowBitOffset = 12;

	// Number of weights
	vertex.positionAndNormal.numWeightsAndBoneIndexes |= ((weights.count - 1) & 0x2) << 30;

	for ( unsigned i = 0; i < weights.count; i++ )
	{
		// bone index
		vertex.positionAndNormal.numWeightsAndBoneIndexes |= (boneIndexes.at (weights.influencers[i]) & 0x1f) << bitOffset;
		bitOffset += 5;

		// overflow of weight
		int weight = static_cast<int>(weights.weights[i] * 1023);

		// Top 2 bits of the weight
		vertex.positionAndNormal.numWeightsAndBoneIndexes |= (weight & 0x300) << overflowBitOffset;
		overflowBitOffset += 2;

		// Bottom 8 bits of the weight
		vertex.positionAndNormal.boneWeightings[i] = static_cast<unsigned char>(weight & 0xff);
	}

	for ( unsigned i = weights.count; i < 4; i++ )
	{
		vertex.positionAndNormal.boneWeightings[i] = 0;
	}
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
	for ( const auto& surface : hierarchy )
	{
		filesize += surface.metadataSize;
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

std::size_t CalculateLODSize ( const ModelDetailData& modelDetail )
{
	std::size_t filesize = 0;

	// LOD data
	const std::vector<GLMSurface>& surfaces = modelDetail.surfaces;

	// Next LOD offset
	filesize += sizeof (int);

	// Surface offsets
	filesize += sizeof (int) * modelDetail.surfaces.size();

	for ( const auto& surface : surfaces )
	{
		filesize += CalculateSurfaceSize (surface.metadata);
	}

	return filesize;
}

std::size_t CalculateLODSize ( const std::vector<ModelDetailData>& modelData )
{
	std::size_t filesize = 0;

	for ( const auto& modelDetail : modelData )
	{
		filesize += CalculateLODSize (modelDetail);
	}

	return filesize;
}

bool GetSurfaceWeightsData (
	const FbxMesh& mesh,
	std::vector<Weights>& vertexWeights,
	std::set<std::string>& referencedBoneNamesSet,
	std::map<std::string, int>& referencedBoneNamesIndex )
{
	for ( auto& weight : vertexWeights )
	{
		weight.count = 0;
	}

	referencedBoneNamesIndex.clear();
	referencedBoneNamesSet.clear();

	FbxSkin *deformer = static_cast<FbxSkin *>(mesh.GetDeformer (0, FbxDeformer::eSkin));
	for ( int j = 0; j < deformer->GetClusterCount(); j++ )
	{
		FbxCluster *cluster = deformer->GetCluster (j);
		int weightCount = cluster->GetControlPointIndicesCount();
		double *weights = cluster->GetControlPointWeights();
		int *weightIndices = cluster->GetControlPointIndices();

		FbxNode *link = cluster->GetLink();
		std::string influencer = link->GetName();

		referencedBoneNamesSet.insert (influencer);

		for ( int k = 0; k < weightCount; k++ )
		{
			Weights& w = vertexWeights[weightIndices[k]];
			w.weights[w.count] = static_cast<float>(weights[k]);
			w.influencers[w.count] = influencer;

			w.count++;
			assert (w.count <= 4);
		}
	}

	int weightCount = 0;

	referencedBoneNamesIndex.clear();
	for ( auto& boneName : referencedBoneNamesSet )
	{
		referencedBoneNamesIndex[boneName] = weightCount++;
	}

	return true;
}

static bool CompareVec3 ( const float v1[3], const FbxVector4& v2, float epsilon )
{
	return fabs (v1[0] - static_cast<float>(v2[0])) < epsilon &&
			fabs (v1[1] - static_cast<float>(v2[1])) < epsilon &&
			fabs (v1[2] - static_cast<float>(v2[2])) < epsilon;
}

static bool CompareVec2 ( const float v1[2], const FbxVector2& v2, float epsilon )
{
	return fabs (v1[0] - static_cast<float>(v2[0])) < epsilon &&
			fabs (v1[1] - static_cast<float>(v2[1])) < epsilon;
}

static int GetNextVertexID (
	const FbxVector4& position,
	const FbxVector4& normal,
	const FbxVector2& texcoord,
	const std::vector<Vertex>& vertices )
{
	for ( size_t i = 0, end = vertices.size(); i < end; i++ )
	{
		if ( CompareVec3 (vertices[i].positionAndNormal.position, position, 1e-3f) &&
			CompareVec3 (vertices[i].positionAndNormal.normal, normal, 1e-6f) &&
			CompareVec2 (vertices[i].texcoord.st, texcoord, 1e-6f) )
		{
			return static_cast<int>(i);
		}
	}

	return static_cast<int>(vertices.size());
}

std::size_t CalculateGLMFileSize (
	const SurfaceHierarchyList& hierarchy,
	const std::vector<ModelDetailData>& modelData )
{
	return CalculateHeaderSize () + 
		CalculateSurfaceHierarchySize (hierarchy) +
		CalculateLODSize (modelData);
}

static bool GetFbxVertexData (
	FbxMesh& mesh,
	int polygon,
	int vertex,
	const FbxVector4 *positions,
	float scale,
	const FbxAMatrix& globalMatrix,
	const char *uvNameSet,
	const std::string& surfaceName,
	FbxVector4& position,
	FbxVector4& normal,
	FbxVector2& texcoord,
	int *positionIndex)
{
	int index = mesh.GetPolygonVertex (polygon, vertex);
	const FbxVector4& fbxPosition = positions[index];
	bool noTexcoord;
	
	*positionIndex = index;
	position = globalMatrix.MultT (fbxPosition) * scale;

	if ( !mesh.GetPolygonVertexNormal (polygon, vertex, normal) )
	{
		std::cerr << "No vertex normal in triangle " << polygon << ", vertex " << vertex << " of " << surfaceName << '\n';
		return false;
	}

	if ( !mesh.GetPolygonVertexUV (polygon, vertex, uvNameSet, texcoord, noTexcoord) || noTexcoord )
	{
		std::cerr << "No texture coordinate in triangle " << polygon << ", vertex " << vertex << " of " << surfaceName << '\n';
		return false;
	}

	return true;
}

ModelDetailData CreateModelLod (
	FbxScene& scene,
	const SurfaceHierarchyList& hierarchy,
	int lod,
	const Skeleton *skeleton,
	std::vector<FbxNode *>& nodes )
{
	// Define these variables outside the loop to prevent repeated
	// allocation and freeing of memory.
	std::vector<Vertex> uniqueVertices;
	std::vector<Weights> vertexWeights;
	std::set<std::string> referencedBoneNamesSet;
	std::map<std::string, int> referencedBoneNamesIndex;

	float scale = skeleton != nullptr ? skeleton->scale : 1.0f;

	ModelDetailData data;
	data.lod = lod;
	data.surfaces.resize (hierarchy.size());

	for ( std::size_t i = 0; i < hierarchy.size(); i++ )
	{
		const mdxmSurfHierarchy_t& surfaceHierarchy = *hierarchy[i].metadata;
		mdxmSurface_t& metadata = data.surfaces[i].metadata;

		FbxNode& node = *nodes[i];
		FbxMesh& mesh = *GetFBXMesh (node);

		int numUVs = mesh.GetElementUVCount();
		int numNormals = mesh.GetElementNormalCount();
		int numPolygons = mesh.GetPolygonCount();
		int numPositions = mesh.GetControlPointsCount();
		int numTriangles = 0;
		const FbxVector4 *positions = mesh.GetControlPoints();

		if ( numNormals == 0 )
		{
			std::cerr << "No normals for surface " << i << " (" << surfaceHierarchy.name << ")\n";
			continue;
		}

		vertexWeights.resize (numPositions);
		GetSurfaceWeightsData (mesh, vertexWeights, referencedBoneNamesSet, referencedBoneNamesIndex);

		std::vector<mdxmTriangle_t>& triangles = data.surfaces[i].triangles;
		triangles.reserve (numPolygons);

		FbxAMatrix globalMatrix = scene.GetEvaluator()->GetNodeGlobalTransform (&node);

		uniqueVertices.clear();

		if ( numUVs == 0 )
		{
			// The only surfaces with no UVs are the tags and the "stupidtriangle" surface.

			// Tags are oriented by taking the last vertex as the origin, and calculating
			// the different edge lengths. The mid-length edge is taken to be the Y axis,
			// and shortest to be X axis.
			const int order[] = {1, 2, 0};
			const FbxGeometryElementNormal *normals = mesh.GetElementNormal();

			if ( numPositions != 3 )
			{
				std::cerr << "Tag " << surfaceHierarchy.name << " should have 3 vertex positions, but has " << numPositions << '\n';
				continue;
			}
			
			// There's going to be 3 vertices for a tag.
			uniqueVertices.resize (3);
			triangles.resize (1);

			for ( int j = 0; j < 3; j++ )
			{
				Vertex& vertex = uniqueVertices[j];
				FbxVector4 positionWS = globalMatrix.MultT (positions[order[j]]) * scale;

				CopyVertexData (vertex,
					positionWS,
					normals->GetDirectArray().GetAt (order[j]),
					vertexWeights[order[j]],
					referencedBoneNamesIndex);

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
			FbxStringList uvNameSetList;
			const char *uvNameSet = nullptr;

			mesh.GetUVSetNames (uvNameSetList);
			if ( uvNameSetList.GetCount() == 0 )
			{
				std::cerr << "ERROR: Surface '" << surfaceHierarchy.name << "' has no UV coordinates.\n";
				continue;
			}

			uvNameSet = uvNameSetList[0];

			for ( int polygon = 0; polygon < numPolygons; polygon++ )
			{
				int vertexIndices[12];
				int numPolygonVertices = mesh.GetPolygonSize (polygon);

				if ( numPolygonVertices > 12 )
				{
					std::cerr << "ERROR: Polygon " << polygon << " in " << surfaceHierarchy.name << " LOD " << lod << " has more than 12 vertices in a polygon.\n";
					break;
				}

				for ( int vertex = 0; vertex < numPolygonVertices; vertex++ )
				{
					int positionIndex;
					FbxVector4 position;
					FbxVector4 normal;
					FbxVector2 texcoord;

					if ( !GetFbxVertexData (
						mesh, polygon, vertex,
						positions, scale, globalMatrix, uvNameSet, surfaceHierarchy.name,
						position, normal, texcoord, &positionIndex) )
					{
						break;
					}

					int vertexID = GetNextVertexID (position, normal, texcoord, uniqueVertices);
					if ( vertexID == uniqueVertices.size() )
					{
						Vertex newVertex;

						CopyVertexData (newVertex,
							position,
							normal,
							vertexWeights[positionIndex],
							referencedBoneNamesIndex);

						newVertex.texcoord.st[0] = static_cast<float>(texcoord[0]);
						newVertex.texcoord.st[1] = static_cast<float>(texcoord[1]);

						uniqueVertices.push_back (newVertex);
					}

					vertexIndices[vertex] = vertexID;
				}

				for ( int i = 1, end = numPolygonVertices - 1; i < end; i++ )
				{
					// Reverse the winding of the triangle vertices.
					mdxmTriangle_t triangle;
					triangle.indexes[0] = vertexIndices[0];
					triangle.indexes[2] = vertexIndices[i];
					triangle.indexes[1] = vertexIndices[i + 1];

					triangles.push_back (triangle);
	
					numTriangles++;
				}
			}
		}

		std::vector<mdxmVertexTexcoord_t>& texcoords = data.surfaces[i].texcoords;
		std::vector<mdxmVertex_t>& vertices = data.surfaces[i].vertices;

		texcoords.resize (uniqueVertices.size());
		vertices.resize (uniqueVertices.size());

		for ( size_t j = 0, numVerts = uniqueVertices.size(); j < numVerts; j++ )
		{
			vertices[j] = uniqueVertices[j].positionAndNormal;
			texcoords[j].st[0] = uniqueVertices[j].texcoord.st[0];
			texcoords[j].st[1] = 1.0f - uniqueVertices[j].texcoord.st[1];
		}

		// Bone references
		std::vector<int>& boneReferences = data.surfaces[i].boneReferences;
		if ( referencedBoneNamesSet.empty() || skeleton == nullptr )
		{
			boneReferences.resize (1);
			boneReferences[0] = 0;
		}
		else
		{
			boneReferences.resize (referencedBoneNamesSet.size());
			for ( const auto& referencedBone : referencedBoneNamesIndex )
			{
				boneReferences[referencedBone.second] = skeleton->boneNamesToIndex.at (referencedBone.first);
			}
		}

		int ofsBoneReferences = 0;
		ofsBoneReferences += sizeof (mdxmSurface_t);
		ofsBoneReferences += sizeof (mdxmTriangle_t) * numTriangles;
		ofsBoneReferences += sizeof (mdxmVertex_t) * uniqueVertices.size();
		ofsBoneReferences += sizeof (mdxmVertexTexcoord_t) * uniqueVertices.size();

		if ( uniqueVertices.size() > 1000 )
		{
			std::cout << "WARNING: '" << surfaceHierarchy.name << "' in LOD " << lod << " has " << uniqueVertices.size() << " vertices - the engine supports up to 1000 vertices per surface.\n";
		}

		if ( numTriangles * 3 > 6000 )
		{
			std::cout << "WARNING: '" << surfaceHierarchy.name << "' in LOD " << lod << " has " << (numTriangles * 3) << " triangles - the engine supports up to 6000 triangles per surface.\n";
		}

		// Fill in the metadata
		metadata.ident = 0;
		metadata.thisSurfaceIndex = i;
		metadata.ofsTriangles = sizeof (mdxmSurface_t);
		metadata.numTriangles = numTriangles;
		metadata.numVerts = uniqueVertices.size();
		metadata.ofsVerts = metadata.ofsTriangles + sizeof (mdxmTriangle_t) * metadata.numTriangles;
		metadata.numBoneReferences = static_cast<int>(boneReferences.size());
		metadata.ofsBoneReferences = ofsBoneReferences;
		metadata.ofsEnd = CalculateSurfaceSize (metadata);
	}

	return data;
}

void CreateSurfaceNameMapping ( const std::vector<FbxNode *>& nodes, std::map<std::string, int>& mapping )
{
	for ( unsigned i = 0, count = nodes.size(); i < count; i++ )
	{
		mapping[nodes[i]->GetName()] = i;
	}
}

bool ReorderSurfacesToLOD0 ( std::vector<FbxNode *>& nodes, const std::map<std::string, int>& nameRemap )
{
	std::vector<FbxNode *> reorderedNodes (nodes.size(), nullptr);

	for ( unsigned i = 0, count = nodes.size(); i < count; i++ )
	{
		const char *nodeName = nodes[i]->GetName();
		std::string surfaceName (nodeName, nodeName + strlen (nodeName) - 2);

		std::map<std::string, int>::const_iterator it = nameRemap.find (surfaceName);
		if ( it == nameRemap.end() )
		{
			std::cerr << "LOD doesn't have surface " << surfaceName << '\n';
			return false;
		}

		reorderedNodes[it->second] = nodes[i];
	}

	nodes = reorderedNodes;

	return true;
}

void LinearizeSurfacesHelper ( FbxNode& root, std::vector<FbxNode *>& nodes )
{
	nodes.push_back (&root);

	for ( int i = 0; i < root.GetChildCount(); i++ )
	{
		LinearizeSurfacesHelper (*root.GetChild (i), nodes);
	}
}

void LinearizeSurfaces ( FbxNode& root, std::vector<FbxNode *>& nodes )
{
	nodes.clear();
	LinearizeSurfacesHelper (root, nodes);
}

std::vector<ModelDetailData> GetModelData (
	FbxScene& scene,
	std::vector<FbxNode *>& modelRoots,
	const Skeleton *skeleton,
	const SurfaceHierarchyList& hierarchy )
{
	std::vector<ModelDetailData> detailData;
	if ( hierarchy.empty() || modelRoots.empty() )
	{
		return detailData;
	}

	std::vector<FbxNode *> nodes;
	int lodCount = modelRoots.size();

	detailData.reserve (lodCount);

	std::map<std::string, int> nameMapping;

	LinearizeSurfaces (*modelRoots[0], nodes);
	CreateSurfaceNameMapping (nodes, nameMapping);

	detailData.push_back (CreateModelLod (scene, hierarchy, 0, skeleton, nodes));

	for ( int i = 1; i < lodCount; i++ )
	{
		LinearizeSurfaces (*modelRoots[i], nodes);

		if ( !ReorderSurfacesToLOD0 (nodes, nameMapping) )
		{
			std::cerr << "Failed to create vertex data for LOD " << i << '\n';

			detailData.push_back (ModelDetailData());

			continue;
		}

		detailData.push_back (CreateModelLod (scene, hierarchy, i, skeleton, nodes));
	}

	return detailData;
}

std::string GetFilename ( const std::string& filePath )
{
	const auto slash = filePath.find_last_of ('/');
	if ( slash == std::string::npos )
	{
		return filePath;
	}

	return filePath.substr (slash + 1);
}

bool WriteLODData ( const mdxmHeader_t& header, char *buffer, const std::vector<ModelDetailData>& modelDetails )
{
	char *lodBasePtr = buffer;

	for ( const auto& modelDetail : modelDetails )
	{
		// Write offset to next header
		int lodOffset = CalculateLODSize (modelDetail);
		int *nextLOD = reinterpret_cast<int *>(lodBasePtr);
		*nextLOD = lodOffset;

		int *surfaceOffsets = reinterpret_cast<int *>(lodBasePtr + sizeof (int));
		char *surfaceData = reinterpret_cast<char *>(surfaceOffsets);
		int surfaceOffset = sizeof (int) * header.numSurfaces;

		for ( int i = 0; i < header.numSurfaces; i++ )
		{
			const GLMSurface& glmSurface = modelDetail.surfaces[i];

			surfaceOffsets[i] = surfaceOffset;

			mdxmSurface_t *surface = reinterpret_cast<mdxmSurface_t *>(&surfaceData[surfaceOffset]);
			*surface = glmSurface.metadata;
			surface->ofsHeader = reinterpret_cast<const char *>(&header) - reinterpret_cast<char *>(surface);
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

	return true;
}

bool WriteSurfaceHierarchyData (
	const mdxmHeader_t& header,
	const SurfaceHierarchyList& surfaceHierarchy,
	int *hierarchyOffsets,
	mdxmSurfHierarchy_t *hierarchy )
{
	char *hierarchyOffsetsBasePtr = reinterpret_cast<char *>(hierarchyOffsets);
	unsigned hierarchyOffset = sizeof (int) * header.numSurfaces;

	for ( int i = 0; i < header.numSurfaces; i++ )
	{
		const GLMSurfaceHierarchy& surface = surfaceHierarchy[i];

		hierarchyOffsets[i] = hierarchyOffset;

		std::memcpy (hierarchy, surface.metadata, surface.metadataSize);

		hierarchyOffset += surface.metadataSize;
		hierarchy = reinterpret_cast<mdxmSurfHierarchy_t *>(&hierarchyOffsetsBasePtr[hierarchyOffset]);
	}

	return true;
}

bool WriteDataToFile ( const std::string& path, const char *begin, const char *end )
{
	std::ofstream file (path.c_str(), std::ios::binary);
	if ( !file )
	{
		return false;
	}

	file.write (begin, end - begin);

	return true;
}

void PrintModelStatistics ( const std::vector<ModelDetailData>& modelDetails )
{
	std::cout << "..." << modelDetails.size() << " LODs\n";
	for ( const auto modelData : modelDetails )
	{
		unsigned numVerts = 0;
		unsigned numTriangles = 0;

		std::cout << "......LOD " << modelData.lod << ": ";
		for ( const auto surface : modelData.surfaces )
		{
			numVerts += surface.metadata.numVerts;
			numTriangles += surface.metadata.numTriangles;
		}

		std::cout << numVerts << " vertices, " << numTriangles << " triangles\n";
	}
}

bool MakeGLMFile (
	FbxScene& scene,
	std::vector<FbxNode *>& modelRoots,
	const Skeleton *skeleton,
	const std::string& outputPath )
{
	if ( modelRoots.empty() )
	{
		return false;
	}

	// Create surface hierarchy
	SurfaceHierarchyList surfaceHierarchy (CreateSurfaceHierarchy (*modelRoots[0]));

	std::cout << "..." << surfaceHierarchy.size() << " surfaces\n";

	// Create surface data
	std::vector<ModelDetailData> modelDetails (GetModelData (scene, modelRoots, skeleton, surfaceHierarchy));

	PrintModelStatistics (modelDetails);

	std::size_t filesize = CalculateGLMFileSize (surfaceHierarchy, modelDetails);

	// Write data
	std::vector<char> buffer (filesize, '\0');

	int hierarchyOffsetsBase = sizeof (mdxmHeader_t);
	int hierarchyBase = hierarchyOffsetsBase + sizeof (int) * surfaceHierarchy.size();
	int lodBase = hierarchyOffsetsBase + CalculateSurfaceHierarchySize (surfaceHierarchy);

	std::string outputFile (GetFilename (outputPath));

	mdxmHeader_t *header = reinterpret_cast<mdxmHeader_t *>(&buffer[0]);
	header->ident = MDXM_IDENT;
	header->version = MDXM_VERSION;
	CopyString (header->name, outputFile.c_str(), sizeof (header->name));
	CopyString (header->animName, ( skeleton == nullptr ) ? "*default" : skeleton->name.c_str(), sizeof (header->animName));
	header->animIndex = -1;
	header->numBones = ( skeleton == nullptr ) ? 1 : skeleton->boneNamesToIndex.size();
	header->numLODs = static_cast<int>(modelDetails.size());
	header->ofsLODs = lodBase;
	header->numSurfaces = static_cast<int>(surfaceHierarchy.size());
	header->ofsSurfHierarchy = hierarchyBase;
	header->ofsEnd = filesize;

	// Surface offsets and hierarchy
	int *hierarchyOffsets = reinterpret_cast<int *>(&buffer[hierarchyOffsetsBase]);
	mdxmSurfHierarchy_t *hierarchy = reinterpret_cast<mdxmSurfHierarchy_t *>(&hierarchyOffsets[header->numSurfaces]);

	if ( !WriteSurfaceHierarchyData (*header, surfaceHierarchy, hierarchyOffsets, hierarchy) )
	{
		std::cerr << "Failed to write surface hierarchy data to file.\n";

		return false;
	}

	// LODs
	if ( !WriteLODData (*header, &buffer[lodBase], modelDetails) )
	{
		std::cerr << "Failed to write vertex data to file.\n";

		return false;
	}

	for ( const auto& surface : surfaceHierarchy )
	{
		char *data = reinterpret_cast<char *>(surface.metadata);
		delete [] data;
	}

	if ( !WriteDataToFile (outputPath, buffer.data(), buffer.data() + buffer.size()) )
	{
		return false;
	}

	return true;
}

void PrintUsage ( const std::string& applicationName )
{
	std::cout << "Usage: " << applicationName << " [-o <GLM file>] <FBX file>\n\n";
}

bool StartsWith ( const std::string& str, const std::string& substr )
{
	return substr.compare (0, substr.length(), str) == 0;
}

std::vector<FbxNode *> GetLodRootModels ( FbxNode& sceneRoot )
{
	// Slightly over size it, but no matter...
	std::vector<FbxNode *> roots;
	roots.reserve (sceneRoot.GetChildCount());

	// Get all the meshes
	for ( int i = 0; i < sceneRoot.GetChildCount(); i++ )
	{
		FbxNode *child = sceneRoot.GetChild (i);

		if ( child->GetNodeAttribute()->GetAttributeType() == FbxNodeAttribute::eMesh )
		{
			roots.push_back (child);
		}
	}

	return roots;
}

} // namespace

int main ( int argc, char *argv[] )
{
	std::vector<std::string> args (argv, argv + argc);
	if ( args.size() == 1 )
	{
		PrintUsage (args[0]);
		return EXIT_FAILURE;
	}

	std::string animationFile;
	std::string outputPath ("model.glm");
	unsigned i;
	for ( i = 1; i < args.size(); i++ )
	{
		if ( args[i].compare ("-o") == 0 )
		{
			i++;
			if ( args.size() > i )
			{
				outputPath = args[i];

				continue;
			}
			else
			{
				PrintUsage (args[0]);

				return EXIT_FAILURE;
			}
		}
		else if ( args[i].compare ("-anim") == 0 )
		{
			i++;
			if ( args.size() > i )
			{
				animationFile = args[i];

				continue;
			}
			else
			{
				PrintUsage (args[0]);

				return EXIT_FAILURE;
			}
		}
	}

	// auto_ptr is deprecated in C++11, but FBX doesn't play nicely with
	// C++11 libraries on OS X, so we have to resort to using auto_ptr
	// instead.
	std::auto_ptr<Skeleton> skeleton;
	std::string modelPath (args.back());

	std::cout << "Converting " << modelPath << " to GLM.\n";

	if ( !animationFile.empty() )
	{
		skeleton = std::auto_ptr<Skeleton>(LoadGLA (animationFile));
		if ( skeleton.get() == nullptr )
		{
			return EXIT_FAILURE;
		}

		std::cout << "Using " << animationFile << " for skeleton.\n";
	}

	FbxManager *fbxManager = FbxManager::Create();
	FbxIOSettings *ios = FbxIOSettings::Create (fbxManager, IOSROOT);
	fbxManager->SetIOSettings (ios);

	FbxImporter *importer = FbxImporter::Create (fbxManager, "");

	if ( !importer->Initialize (modelPath.c_str(), -1, fbxManager->GetIOSettings()) )
	{
		std::cerr << "Failed to import " << modelPath << ".\n";

		fbxManager->Destroy();

		return EXIT_FAILURE;
	}

	FbxScene *scene = FbxScene::Create (fbxManager, "model");
	importer->Import (scene);
	importer->Destroy();

	FbxNode *root = scene->GetRootNode();
	if ( root == nullptr )
	{
		std::cerr << "No root node :(\n";

		fbxManager->Destroy();

		return EXIT_FAILURE;
	}

	std::vector<FbxNode *> modelRoots (GetLodRootModels (*root));
	std::sort (modelRoots.begin(), modelRoots.end(),
		[]( const FbxNode *a, const FbxNode *b ) { return std::strcmp (a->GetName(), b->GetName()) < 0; });

	if ( !MakeGLMFile (*scene, modelRoots, skeleton.get(), outputPath) )
	{
		std::cerr << "Failed to create GLM file " << outputPath << ".\n";

		fbxManager->Destroy();

		return EXIT_FAILURE;
	}

	std::cout << "GLM file has been written to " << outputPath << ".\n";

	fbxManager->Destroy();

	return 0;
}
