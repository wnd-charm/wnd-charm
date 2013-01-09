/* 
 * hashlib++ - a simple hash library for C++
 * 
 * Copyright (c) 2007-2010 Benjamin Grüdelbach
 * 
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * 	1)     Redistributions of source code must retain the above copyright
 * 	       notice, this list of conditions and the following disclaimer.
 * 
 * 	2)     Redistributions in binary form must reproduce the above copyright
 * 	       notice, this list of conditions and the following disclaimer in
 * 	       the documentation and/or other materials provided with the
 * 	       distribution.
 * 	     
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

//----------------------------------------------------------------------	

/*
 * The hashlib++ MD5 implementation is derivative from the sourcecode
 * published in RFC 1321 
 * 
 * Copyright (C) 1991-2, RSA Data Security, Inc. Created 1991. All
 * rights reserved.
 * 
 * License to copy and use this software is granted provided that it
 * is identified as the "RSA Data Security, Inc. MD5 Message-Digest
 * Algorithm" in all material mentioning or referencing this software
 * or this function.
 * 
 * License is also granted to make and use derivative works provided
 * that such works are identified as "derived from the RSA Data
 * Security, Inc. MD5 Message-Digest Algorithm" in all material
 * mentioning or referencing the derived work.
 * 
 * RSA Data Security, Inc. makes no representations concerning either
 * the merchantability of this software or the suitability of this
 * software for any particular purpose. It is provided "as is"
 * without express or implied warranty of any kind.
 * 
 * These notices must be retained in any copies of any part of this
 * documentation and/or software.
 */

//----------------------------------------------------------------------	

/**
 *  @file 	md5.h
 *  @brief	This file contains the declaration of the MD5 class
 *  @date 	Mo 17 Sep 2007
 */  

//---------------------------------------------------------------------- 
//include protection
#ifndef MD5_H
#define MD5_H

//---------------------------------------------------------------------- 
//hl includes
// use standard types in <stdint.h>
#include <stdint.h>
#include <string>

namespace MD5 {

//---------------------------------------------------------------------- 
//typedefs
typedef uint8_t *POINTER;

/**
 * @brief this struct represents a MD5-hash context.
 */
typedef struct 
{
	/** state (ABCD) */
	unsigned long int state[4];   	      

	/** number of bits, modulo 2^64 (lsb first) */
	unsigned long int count[2];

	/** input buffer */
	unsigned char buffer[64];
} HL_MD5_CTX;

//---------------------------------------------------------------------- 
//defines
const static size_t HashSize=16;

//---------------------------------------------------------------------- 

/**
 *  @brief 	This class represents the implementation of 
 *   		the md5 message digest algorithm.
 *
 *   		Basically the class provides three public member-functions
 *   		to create a hash:  Init(), Update() and Final().
 */  
class MD5
{

	private:

		/**
		 *  @brief 	Basic transformation. Transforms state based on block.
		 *  @param	state	state to transform
		 *  @param	block	block to transform
		 */  
		void Transform (unsigned long int state[4], const unsigned char block[64]);

		/**
		 *  @brief 	Encodes input data
		 *  @param	output Encoded data as OUT parameter
		 *  @param	input Input data
		 *  @param	len The length of the input assuming it is a
		 *  		multiple of 4
		 */  
		void Encode (unsigned char* output,
			     const unsigned long int *input,
			     size_t len);

		/**
		 *  @brief 	Decodes input data into output
		 *  @param	output Decoded data as OUT parameter
		 *  @param	input Input data
		 *  @param	len The length of the input assuming it is a
		 *  		multiple of 4
		 */  
		void Decode (unsigned long int *output,
			     const unsigned char *input,
			     size_t len);

			HL_MD5_CTX context;
	public:
	
		/**
		 *  @brief 	constructor - calls Reset
		 */  
		MD5( ) {Reset();};

		/**
		 *  @brief 	Initialization begins an operation,
		 *  		writing a new context
		 */  
		void Reset ();

		/**
		 *  @brief 	Block update operation. Continues an md5
		 *  		message-digest operation, processing another
		 *  		message block, and updating the context.
		 *  @param	input The data to write into the context
		 *  @param	inputLen The length of the input data
		 */  
		void Update (const unsigned char *input,
			       	const size_t inputLen);

		/**
		 *  @brief 	Finalization ends the md5 message-digest 
		 *  		operation, writing the the message digest and
		 *  		zeroizing the context.
		 *  @param	digest This is an OUT parameter which contains
		 *  		the created hash after the method returns
		 */  
		void Result (unsigned char digest[HashSize]);
};



/**
 *  @brief 	This computes a md5 digest on a string in one go.
 *
 *  @param	input_str the std::string input string
 *  @return	an std::string of HashSize bytes (16)
 */  
const std::string &fromString( const std::string &input_str, std::string &output_str) ;

/**
 *  @brief 	This computes a md5 digest on a set of bytes in one go.
 *
 *  @param	bytes pointer to bytes
 *  @param	length number of bytes pointed to by the bytes pointer
 *  @return	an std::string of HashSize bytes (16)
 */  
const std::string &fromBytes( const unsigned char *bytes, const size_t length, std::string &output_str) ;

//----------------------------------------------------------------------
//end of MD5 namespace
}

//---------------------------------------------------------------------- 
//End of include protection
#endif

//---------------------------------------------------------------------- 
//EOF
