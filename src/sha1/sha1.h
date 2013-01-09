/* 
 * hashlib++ - a simple hash library for C++
 * 
 * Copyright (c) 2007-2010 Benjamin Gr�delbach
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
 * The hashlib++ SHA1 implementation is derivative from the sourcecode
 * published in RFC 3174  
 *
 * Copyright (C) The Internet Society (2001).  All Rights Reserved.
 * 
 * This document and translations of it may be copied and furnished to
 * others, and derivative works that comment on or otherwise explain it
 * or assist in its implementation may be prepared, copied, published
 * and distributed, in whole or in part, without restriction of any
 * kind, provided that the above copyright notice and this paragraph are
 * included on all such copies and derivative works.  However, this
 * document itself may not be modified in any way, such as by removing
 * the copyright notice or references to the Internet Society or other
 * Internet organizations, except as needed for the purpose of
 * developing Internet standards in which case the procedures for
 * copyrights defined in the Internet Standards process must be
 * followed, or as required to translate it into languages other than
 * English.
 * 
 * The limited permissions granted above are perpetual and will not be
 * revoked by the Internet Society or its successors or assigns.
 * 
 * This document and the information contained herein is provided on an
 * "AS IS" basis and THE INTERNET SOCIETY AND THE INTERNET ENGINEERING
 * TASK FORCE DISCLAIMS ALL WARRANTIES, EXPRESS OR IMPLIED, INCLUDING
 * BUT NOT LIMITED TO ANY WARRANTY THAT THE USE OF THE INFORMATION
 * HEREIN WILL NOT INFRINGE ANY RIGHTS OR ANY IMPLIED WARRANTIES OF
 * MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.
 */

//---------------------------------------------------------------------- 

/**
 *  @file 	hl_sha1.h
 *  @brief	This file contains the declaration of the SHA1 class
 *  @date 	Mo 17 Sep 2007
 */  

//---------------------------------------------------------------------- 
//include protection
#ifndef SHA1_H
#define SHA1_H

//---------------------------------------------------------------------- 
//hl includes
// use standard types in <stdint.h>
#include <stdint.h>
#include <string>

namespace SHA1 {
//---------------------------------------------------------------------- 
//enums

#ifndef _SHA_enum_
#define _SHA_enum_
enum
{
    shaSuccess = 0,
    shaNull,            /* Null pointer parameter */
    shaInputTooLong,    /* input data too long */
    shaStateError       /* called Input after Result */
};
#endif

//---------------------------------------------------------------------- 
//defines
const static size_t HashSize=20;

//---------------------------------------------------------------------- 
//structs

/**
 * @brief this struct represents a SHA1-hash context.
 */
typedef struct HL_SHA1_CTX
{
	/** Message Digest */
	uint32_t Intermediate_Hash[HashSize/4];

	/** Message length in bits */
	uint32_t Length_Low;            

	/** Message length in bits */
	uint32_t Length_High;

	/** Index into message block array */
	uint16_t Message_Block_Index;

	/** 512-bit message blocks */
	uint8_t Message_Block[64];      

	/** Is the digest computed? */
	int Computed;

	/** Is the message digest corrupted? */
	int Corrupted;

} HL_SHA1_CTX;

//---------------------------------------------------------------------- 
//class definition

/**
 *  @brief 	This class represents the implementation of 
 *   		the sha1 algorithm.
 *
 *   		Basically the class provides three public member-functions
 *   		to create a hash:  SHA1Reset(), SHA1Input() and SHA1Result().
 *   		If you want to create a hash based on a string or file quickly
 *   		you should use the sha1wrapper class instead of SHA1.
 */  
class SHA1
{
	private:

			/**
			 *  @brief 	Internal method to padd the message
			 *
			 *      	According to the standard, the message must
			 *      	be padded to an even 512 bits. The first 
			 *      	padding bit must be a '1'.  The last 64	bits 
			 *      	represent the length of the original message.
			 *      	All bits in between should be 0.
			 *      	This function will pad the message according 
			 *      	to those rules by filling the Message_Block array
			 *      	accordingly.  It will also call the 
			 *      	ProcessMessageBlock function provided appropriately.
			 *      	When it returns, it can be assumed that the message
			 *      	digest has been computed.
			 *
			 *
			 */  
			void PadMessage();

			/**
			 *  @brief      This member-function will process the next 512 bits of the
			 *  		message stored in the Message_Block array.
			 *
			 *      	Many of the variable names in this code, especially the
			 *      	single character names, were used because those were the
			 *      	names used in the publication.
			 *
			 */  
			void ProcessMessageBlock();
			HL_SHA1_CTX context;

	public:

		/**
		 *  @brief 	constructor - calls Reset
		 */  
		SHA1( ) {Reset();};


		/**
		 *  @brief 	Resets the sha1 context and starts a new
		 *  		hashprocess
		 *  @return	0 on succes an error number otherwise
		 */  
		int Reset( );

		/**
		 *  @brief 	Data input.
		 *
		 *  		This memberfunction add data to the specified
		 *  		context.
		 *
		 *  @param	message_array The data to add
		 *  @param	length The length of the data to add
		 */  
		int Input(  const uint8_t *message_array,
				size_t  length);

		/**
		 *  @brief 	This ends the sha operation, zeroizing the context
		 *  		and returning the computed hash.
		 *
		 *  @param	Message_Digest This is an OUT parameter which
		 *  		contains the hash after the menberfunction returns
		 *  @return	0 on succes, an error-code otherwise
		 */  
		int Result( uint8_t     Message_Digest[HashSize]);
};



/**
 *  @brief 	This computes a sha1 digest on a string in one go.
 *
 *  @param	input_str the std::string input string
 *  @return	an std::string of SHA1HashSize bytes (20)
 */  
const std::string &fromString( const std::string &input_str, std::string &output_str) ;

/**
 *  @brief 	This computes a sha1 digest on a set of bytes in one go.
 *
 *  @param	bytes pointer to bytes
 *  @param	length number of bytes pointed to by the bytes pointer
 *  @return	an std::string of SHA1HashSize bytes (20)
 */  
const std::string &fromBytes( const unsigned char *bytes, const size_t length, std::string &output_str) ;

//----------------------------------------------------------------------
//end of SHA1 namespace
}

//----------------------------------------------------------------------
//end of include protection
#endif

//----------------------------------------------------------------------
//EOF
