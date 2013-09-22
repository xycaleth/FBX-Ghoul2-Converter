#include <algorithm>
#include <cassert>
#include <cctype>
#include <fbxsdk/fileio/fbxiosettings.h>
#include <fstream>
#include <iostream>
#include <fbxsdk.h>
#include <map>
#include <vector>

#include "GLM.h"

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

FbxMesh *GetFBXMesh ( FbxNode& node, int lod )
{
	// Default attribute first
	if ( lod == 0 )
	{
		FbxNodeAttribute *defaultAttribute = node.GetNodeAttribute();
		if ( defaultAttribute != nullptr )
		{
			return GetFBXMesh (*defaultAttribute);
		}

		return nullptr;
	}

	FbxMesh *mesh = nullptr;
	for ( int i = 0; i < node.GetNodeAttributeCount(); i++ )
	{
		FbxNodeAttribute *attribute = node.GetNodeAttributeByIndex (i);
		FbxMesh *lodMesh = GetFBXMesh (*attribute);
		if ( lodMesh != nullptr )
		{
			mesh = lodMesh;
			lod--;

			if ( lod == 0 )
			{
				break;
			}
		}
	}

	if ( lod > 0 )
	{
		return nullptr;
	}

	return mesh;
}

FbxNode *GetRootMesh ( FbxNode& root )
{
	for ( int i = 0, count = root.GetChildCount(); i < count; i++ )
	{
		FbxNode *child = root.GetChild (i);
		FbxNodeAttribute *defaultAttribute = child->GetNodeAttribute();
		if ( defaultAttribute == nullptr )
		{
			continue;
		}

		if ( defaultAttribute->GetAttributeType() == FbxNodeAttribute::eMesh )
		{
			return child;
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
		CopyString (&hierarchyNode->name[1], node.GetName() + 5, sizeof (hierarchyNode->name) - 1);
	}
	else
	{
		CopyString (hierarchyNode->name, node.GetName() + 5, sizeof (hierarchyNode->name));
	}
	std::transform (hierarchyNode->name, hierarchyNode->name + sizeof (hierarchyNode->name),
						hierarchyNode->name, []( char c ) { return std::tolower (c); });

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

int GetLodCount ( FbxNode& rootMesh )
{
	int lod = 0;
	while ( GetFBXMesh (rootMesh, lod) != nullptr )
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

void CopyVertexData ( Vertex& vertex, const FbxAMatrix& globalMatrix, const FbxVector4& position, const FbxVector4& normal )
{
	FbxVector4 positionWS = globalMatrix.MultT (position);
	vertex.positionAndNormal.position[0] = static_cast<float>(positionWS[0]);
	vertex.positionAndNormal.position[1] = static_cast<float>(positionWS[1]);
	vertex.positionAndNormal.position[2] = static_cast<float>(positionWS[2]);

	FbxVector4 normalWS = globalMatrix.MultT (normal);
	vertex.positionAndNormal.normal[0] = static_cast<float>(normalWS[0]);
	vertex.positionAndNormal.normal[1] = static_cast<float>(normalWS[1]);
	vertex.positionAndNormal.normal[2] = static_cast<float>(normalWS[2]);

	vertex.positionAndNormal.numWeightsAndBoneIndexes = 0u;
	vertex.positionAndNormal.boneWeightings[0] = 255u;
	vertex.positionAndNormal.boneWeightings[1] = 0u;
	vertex.positionAndNormal.boneWeightings[2] = 0u;
	vertex.positionAndNormal.boneWeightings[3] = 0u;
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

std::size_t CalculateLODSize (
		const std::vector<ModelDetailData>& modelData )
{
	std::size_t filesize = 0;

	for ( const auto& modelDetail : modelData )
	{
		filesize += CalculateLODSize (modelDetail);
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

ModelDetailData CreateModelLod ( FbxScene& scene, const SurfaceHierarchyList& hierarchy, int lod )
{
	ModelDetailData data;
	data.lod = lod;
	data.surfaces.resize (hierarchy.size());

	// Define these variables outside the loop to prevent repeated
	// allocation and freeing of memory.
	std::map<VertexId, int> uniqueVerticesMap;
	std::vector<Vertex> uniqueVertices;

	for ( std::size_t i = 0; i < hierarchy.size(); i++ )
	{
		const mdxmSurfHierarchy_t& surfaceHierarchy = *hierarchy[i].metadata;

		FbxNode& fbxNode = *hierarchy[i].node;
		mdxmSurface_t& metadata = data.surfaces[i].metadata;

		FbxMesh& mesh = *GetFBXMesh (fbxNode, lod);
		const int *indices = mesh.GetPolygonVertices();
		const FbxVector4 *positions = mesh.GetControlPoints();
		const FbxLayer *layer0 = mesh.GetLayer (0);

		if ( layer0 == nullptr )
		{
			std::cerr << "Normals and texture coordinates not found for surface " << i << " (" << surfaceHierarchy.name << ")\n";
		}

		const FbxLayerElementNormal *normalsLayer = layer0->GetNormals();
		if ( normalsLayer == nullptr )
		{
			std::cerr << "No normals for surface " << i << " (" << surfaceHierarchy.name << ")\n";
		}

		FbxLayerElementArrayTemplate<FbxVector4>& normals = normalsLayer->GetDirectArray();

		const FbxLayerElementUV *uvLayer = layer0->GetUVs();
		FbxLayerElementArrayTemplate<FbxVector2> *uvs = nullptr;

		std::size_t numUniquePoints = mesh.GetControlPointsCount();
		std::size_t numPositions = mesh.GetControlPointsCount();
		if ( uvLayer != nullptr )
		{
			if ( !mesh.GetTextureUV (&uvs) )
			{
				std::cerr << "Failed to retrieve texture UVs for surface " << i << '(' << surfaceHierarchy.name << ")\n";
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
			uniqueVertices.reserve (numUniquePoints);
		}

		FbxAMatrix globalMatrix;
		globalMatrix = scene.GetEvaluator()->GetNodeGlobalTransform (&fbxNode);

		if ( uvs == nullptr )
		{
			// The only surfaces with no UVs are the tags and the "stupidtriangle" surface.
			
			// Tags are oriented by taking the last vertex as the origin, and calculating
			// the different edge lengths. The mid-length edge is taken to be the Y axis,
			// and shortest to be X axis.
			const int order[] = {1, 2, 0};
			for ( int j = 0; j < 3; j++ )
			{
				Vertex vertex;

				CopyVertexData (vertex, globalMatrix, positions[order[j]], normals[order[j]]);

				// Not strictly necessary but makes things cleaner in the output file.
				vertex.texcoord.st[0] = 0.0f;
				vertex.texcoord.st[1] = 0.0f;

				uniqueVertices.push_back (vertex);
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
				// With the same number of texture coordinates, and the same number of
				// position/normals, it's much easier to produce unique vertices.
				uniqueVertices.resize (numUniquePoints);
				for ( int j = 0, count = uvIndices.GetCount(); j < count; j++ )
				{
					Vertex& vertex = uniqueVertices[indices[j]];
					if ( vertex.filled )
					{
						continue;
					}

					vertex.filled = true;

					CopyVertexData (vertex, globalMatrix, positions[indices[j]], normals[indices[j]]);

					FbxVector2 tc = uvs->GetAt (uvIndices[j]);
					vertex.texcoord.st[0] = static_cast<float>(tc[0]);
					vertex.texcoord.st[1] = 1.0f - static_cast<float>(tc[1]);
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
				int index = 0;

				uniqueVerticesMap.clear();
				for ( int tri = 0, k = 0; tri < numTriangles; tri++ )
				{
					int triangle[3];

					for ( int j = 0; j < 3; k++, j++ )
					{
						VertexId id (indices[k], uvIndices[k]);
						const auto it = uniqueVerticesMap.find (id);

						if ( it == uniqueVerticesMap.end() )
						{
							Vertex v;

							uniqueVerticesMap.insert (std::make_pair (id, index));
							CopyVertexData (v, globalMatrix, positions[id.positionId], normals[id.positionId]);

							FbxVector2 tc = uvs->GetAt (id.texcoordId);
							v.texcoord.st[0] = static_cast<float>(tc[0]);
							v.texcoord.st[1] = 1.0f - static_cast<float>(tc[1]);

							triangle[j] = index;
							index++;

							uniqueVertices.push_back (v);
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
		std::vector<mdxmVertex_t>& vertices = data.surfaces[i].vertices;

		texcoords.resize (numVerts);
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

		int ofsBoneReferences = 0;
		ofsBoneReferences += sizeof (mdxmSurface_t);
		ofsBoneReferences += sizeof (mdxmTriangle_t) * metadata.numTriangles;
		ofsBoneReferences += sizeof (mdxmVertex_t) * metadata.numVerts;
		ofsBoneReferences += sizeof (mdxmVertexTexcoord_t) * metadata.numVerts;

		// Fill in the metadata
		metadata.ident = 0;
		metadata.thisSurfaceIndex = i;
		metadata.ofsTriangles = sizeof (mdxmSurface_t);
		metadata.numTriangles = numTriangles;
		metadata.numVerts = numVerts;
		metadata.ofsVerts = metadata.ofsTriangles + sizeof (mdxmTriangle_t) * metadata.numTriangles;
		metadata.numBoneReferences = 1;
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
	for ( int i = 0; i < lodCount; i++ )
	{
		detailData.push_back (CreateModelLod (scene, hierarchy, i));
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

bool WriteLoDData ( const mdxmHeader_t& header, char *buffer, const std::vector<ModelDetailData>& modelDetails )
{
	char *lodBasePtr = buffer;
	int lodOffset = 0;

	for ( const auto& modelDetail : modelDetails )
	{
		// Write offset to next header
		int lodOffset = CalculateLODSize (modelDetail);
		int *nextLOD = reinterpret_cast<int *>(lodBasePtr);
		*nextLOD = lodOffset;

		int *surfaceOffsets = reinterpret_cast<int *>(lodBasePtr + sizeof (int));
		char *surfaceData = reinterpret_cast<char *>(surfaceOffsets);
		int surfaceOffset = sizeof (int) * header.numSurfaces;

		for ( unsigned i = 0; i < header.numSurfaces; i++ )
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

bool WriteSurfaceHierarchyData ( const mdxmHeader_t& header,
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
	std::cout << "..." << modelDetails.size() << " LoDs\n";
	for ( const auto modelData : modelDetails )
	{
		unsigned numVerts = 0u;
		unsigned numTriangles = 0u;

		std::cout << "......LoD " << modelData.lod << ": ";
		for ( const auto surface : modelData.surfaces )
		{
			numVerts += surface.metadata.numVerts;
			numTriangles += surface.metadata.numTriangles;
		}

		std::cout << numVerts << " vertices, " << numTriangles << " triangles\n";
	}
}

bool MakeGLMFile ( FbxScene& scene, FbxNode& root, const std::string& outputPath )
{
	FbxNode *meshRoot = GetRootMesh (root);

	// Create surface hierarchy
	SurfaceHierarchyList surfaceHierarchy (CreateSurfaceHierarchy (*meshRoot));

	std::cout << "..." << surfaceHierarchy.size() << " surfaces\n";

	// Create surface data
	std::vector<ModelDetailData> modelDetails (GetModelData (scene, surfaceHierarchy));

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
	CopyString (header->animName, "*default", sizeof (header->animName));
	header->animIndex = -1;
	header->numBones = 1;
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
	if ( !WriteLoDData (*header, &buffer[lodBase], modelDetails) )
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

int main ( int argc, char *argv[] )
{
	std::vector<std::string> args (argv, argv + argc);
	if ( args.size()  == 1 )
	{
		PrintUsage (args[0]);
		return 1;
	}

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

				return 1;
			}
		}
		else
		{
			break;
		}
	}

	if ( i == args.size() )
	{
		PrintUsage (args[0]);

		return 1;
	}

	std::string modelPath (args[i]);
	std::cout << "Converting " << modelPath << " to GLM.\n";

	FbxManager *fbxManager = FbxManager::Create();
	FbxIOSettings *ios = FbxIOSettings::Create (fbxManager, IOSROOT);
	fbxManager->SetIOSettings (ios);

	FbxImporter *importer = FbxImporter::Create (fbxManager, "");

	if ( !importer->Initialize (modelPath.c_str(), -1, fbxManager->GetIOSettings()) )
	{
		importer->Destroy();
		std::cerr << "Failed to import " << modelPath << ".\n";

		return 1;
	}

	FbxScene *scene = FbxScene::Create (fbxManager, "model");
	importer->Import (scene);
	importer->Destroy();

	FbxNode *root = scene->GetRootNode();
	if ( root == nullptr )
	{
		std::cerr << "No root node :(\n";

		return 1;
	}

	if ( !MakeGLMFile (*scene, *root, outputPath) )
	{
		std::cerr << "Failed to create GLM file " << outputPath << ".\n";

		return 1;
	}

	std::cout << "GLM file has been written to " << outputPath << ".\n";

	return 0;
}
