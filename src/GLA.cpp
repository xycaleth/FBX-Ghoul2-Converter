#include "GLA.h"

#include <iostream>
#include <vector>

#include "IO.h"

Skeleton *LoadGLA ( const std::string& animationPath )
{
	std::vector<char> buffer;
	if ( !ReadFile (animationPath, buffer) )
	{
		std::cerr << "ERROR: " << animationPath << " could not be read.\n";
		return nullptr;
	}

	mdxaHeader_t *header = reinterpret_cast<mdxaHeader_t *>(&buffer[0]);
	if ( header->ident != MDXA_IDENT )
	{
		std::cerr << "ERROR: GLA header is incorrect (expected " << MDXA_IDENT << ", found " << header->ident << ".\n";
		return nullptr;
	}

	if ( header->version != MDXA_VERSION )
	{
		std::cerr << "ERROR: GLA file has wrong version (expected " << MDXA_VERSION << ", found " << header->version << ".\n";
		return nullptr;
	}

	char *skelBase = &buffer[sizeof (mdxaHeader_t)];
	mdxaSkelOffsets_t *offsets = reinterpret_cast<mdxaSkelOffsets_t *>(skelBase);
	Skeleton *skeleton = new Skeleton();

	skeleton->name = header->name;
	skeleton->scale = header->scale;

	for ( int i = 0; i < header->numBones; i++ )
	{
		mdxaSkel_t *skel = reinterpret_cast<mdxaSkel_t *>(skelBase + offsets->offsets[i]);
		skeleton->boneNamesToIndex.insert (std::make_pair (skel->name, i));
	}

	return skeleton;
}

