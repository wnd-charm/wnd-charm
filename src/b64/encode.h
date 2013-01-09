// :mode=c++:
/*
encode.h - c++ wrapper for a base64 encoding algorithm

This is part of the libb64 project, and has been placed in the public domain.
For details, see http://sourceforge.net/projects/libb64
*/
#ifndef BASE64_ENCODE_H
#define BASE64_ENCODE_H

#include <iostream>

namespace base64
{
	extern "C" 
	{
		#include "cencode.h"
	}
	
	const static int BUFFERSIZE = 64;

	struct encoder
	{
		base64_encodestate _state;
		int _buffersize;

		encoder(int buffersize_in = BUFFERSIZE)
		: _buffersize(buffersize_in)
		{}

		inline size_t calc_encoded_size (const size_t plainlength) {
			size_t coded_size = ((plainlength + ( (plainlength % 3) ? (3 - (plainlength % 3)) : 0) ) / 3) * 4;
			if (_state.chars_per_line)
				coded_size += (((coded_size) / _state.chars_per_line) * 2) + 1;
			return (coded_size);
		}
		int encode(char value_in)
		{
			return base64_encode_value(value_in);
		}

		int encode(const char* code_in, const int length_in, char* plaintext_out)
		{
			return base64_encode_block(code_in, length_in, plaintext_out, &_state);
		}

		int encode_end(char* plaintext_out)
		{
			return base64_encode_blockend(plaintext_out, &_state);
		}

		std::ostream &encode(std::istream& istream_in, std::ostream& ostream_in, const bool newlines = true)
		{
			base64_init_encodestate(&_state);
			if (!newlines) _state.chars_per_line = 0;
			//
			const int N = _buffersize;
			char* plaintext = new char[N];
			char* code = new char[calc_encoded_size (N)];
			size_t plainlength;
			size_t codelength;

			do
			{
				istream_in.read(plaintext, N);
				plainlength = istream_in.gcount();
				//
				codelength = encode(plaintext, plainlength, code);
				ostream_in.write(code, codelength);
			}
			while (istream_in.good() && plainlength > 0);

			codelength = encode_end(code);
			ostream_in.write(code, codelength);
			//
			base64_init_encodestate(&_state);

			delete [] code;
			delete [] plaintext;
			return (ostream_in);
		}

		std::string &encode(const char *plaintext, const size_t plainlength, std::string &out, const bool newlines = true)
		{
			base64_init_encodestate(&_state);
			if (!newlines) _state.chars_per_line = 0;
			//
			int codelength = calc_encoded_size (plainlength);

			// writing directly into std::string.data() is not recommended. Oh well.
			out.resize (codelength);
			char* code = (char *)out.data();
			int codelength_out,codelength_end;

			codelength_out = encode(plaintext, plainlength, code);
			codelength_end = encode_end(code + codelength_out);			
			return (out);
		}

		std::string &encode(const std::string &plaintext, std::string &out, const bool newlines = true)
		{
			base64_init_encodestate(&_state);
			if (!newlines) _state.chars_per_line = 0;
			size_t plainlength = plaintext.length();
			//
			int codelength = calc_encoded_size (plainlength);

			// writing directly into std::string.data() is not recommended. Oh well.
			out.resize (codelength);
			char* code = (char *)out.data();
			int codelength_out,codelength_end;

			codelength_out = encode(plaintext.data(), plainlength, code);
			codelength_end = encode_end(code + codelength_out);			
			return (out);
		}
	};

} // namespace base64

#endif // BASE64_ENCODE_H

