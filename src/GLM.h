#ifndef GLM_H
#define GLM_H

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

const int G2SURFACE_BOLT	= 0x1;
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

	/*!
	 * Packed 32-bit uint. This stores three things: the number of weights,
	 * the indexes of the bones which this vertex is weighted to, and the top
	 * 2 bits of each of the weights influences.
	 *
	 * The uint is laid like like so:
	 * - Bits 31 - 30:	Number of weights (where 0 = 1, and up to 4 weights).
	 * - Bits 29 - 28:	Reserved
	 * - Bits 27 - 26:	Upper 2 bits of weightings of bone 4.
	 * - Bits 25 - 24:	Upper 2 bits of weightings of bone 3.
	 * - Bits 23 - 22:	Upper 2 bits of weightings of bone 2.
	 * - Bits 21 - 20:	Upper 2 bits of weightings of bone 1.
	 * - Bits 19 - 15:	Bone index 4
	 * - Bits 14 - 10:	Bone index 3
	 * - Bits  9 - 5:	Bone index 2
	 * - Bits  4 - 0:	Bone index 1
	 */
	unsigned int numWeightsAndBoneIndexes;

	/*!
	 * Stores the lower 8 bits of the bone weightings. When combined with top 2
	 * bits from the above field, the result should be divided by 1023 to get
	 * the floating-point weight.
	 */
	unsigned char boneWeightings[4];
};

struct mdxmVertexTexcoord_t
{
	float st[2];
};

struct GLMSurfaceHierarchy
{
	mdxmSurfHierarchy_t *metadata;
	std::size_t metadataSize;
};
typedef std::vector<GLMSurfaceHierarchy> SurfaceHierarchyList;

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

#endif
