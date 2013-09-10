#ifndef __SPLIT_H
#define __SPLIT_H

// functions to split a string by a specific delimiter
#include <string>
#include <vector>
#include <sstream>
#include <string.h>

// thanks to Evan Teran, http://stackoverflow.com/questions/236129/how-to-split-a-string/236803#236803

// split a string on a single delimiter character (delim)
std::vector<std::string>& split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string>  split(const std::string &s, char delim);

// split a string on any character found in the string of delimiters (delims)
std::vector<std::string>& split(const std::string &s, const std::string& delims, std::vector<std::string> &elems);
std::vector<std::string>  split(const std::string &s, const std::string& delims);

// from Marius, http://stackoverflow.com/a/1493195/238609
template < class ContainerT >
void tokenize(const std::string& str, ContainerT& tokens,
              const std::string& delimiters = " ", const bool trimEmpty = false)
{

    std::string::size_type pos, lastPos = 0;
    while(true)
    {
	pos = str.find_first_of(delimiters, lastPos);
	if(pos == std::string::npos)
	{

	    pos = str.length();

	    if(pos != lastPos || !trimEmpty) {
		tokens.push_back(typename ContainerT::value_type(str.data()+lastPos, (typename ContainerT::value_type::size_type)pos-lastPos));
	    }

	    break;
	}
	else
	{
	    if(pos != lastPos || !trimEmpty) {
		tokens.push_back(typename ContainerT::value_type(str.data()+lastPos, (typename ContainerT::value_type::size_type)pos-lastPos));
	    }
	}

	lastPos = pos + 1;
    }
};


#endif
