#include "IO.h"

#include <fstream>

bool ReadFile ( const std::string& animationPath, std::vector<char>& buffer )
{
	std::ifstream ifs (animationPath.c_str(), std::ios::binary);
	if ( !ifs )
	{
		return false;
	}

	std::streamoff fileLen;

	ifs.seekg (0, std::ios::end);
	fileLen = ifs.tellg();
	ifs.seekg (0, std::ios::beg);

	buffer.resize (static_cast<unsigned int>(fileLen));
	ifs.read (buffer.data(), fileLen);

	return true;
}
