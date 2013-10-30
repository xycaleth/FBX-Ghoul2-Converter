#ifndef GLA_H
#define GLA_H

const int MDXA_IDENT = ('A' << 24) | ('G' << 16) | ('L' << 8) | '2';
const int MDXA_VERSION = 6;

struct mdxaHeader_t
{
	int ident;
	int version;
	char name[64];
	float scale;
	int numFrames;
	int ofsFrames;
	int numBones;
	int ofsCompBonePool;
	int ofsSkel;
	int ofsEnd;
};

struct mdxaSkelOffsets_t
{
	int offsets[1];
};

struct mdxaBone_t
{
	float matrix[3][4];
};

struct mdxaSkel_t
{
	char name[64];
	unsigned int flags;
	int parent;
	mdxaBone_t basePose;
	mdxaBone_t inverseBasePose;
	int numChildren;
	int children[1];
};

struct Skeleton
{
	std::string name;
	std::map<std::string, int> boneNamesToIndex;
	float scale;
};

#endif
