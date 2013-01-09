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
 *  @file 	sha1.cpp
 *  @brief	This file contains the implementation of the SHA1 class
 *  @date 	Mo 17 Sep 2007
 */  

//---------------------------------------------------------------------- 
//hashlib++ includes
#include "sha1.h"

namespace SHA1 {


//---------------------------------------------------------------------- 
//defines

/*
 *  Define the SHA1 circular left shift macro
 */
#define SHA1CircularShift(bits,word) \
                (((word) << (bits)) | ((word) >> (32-(bits))))

//----------------------------------------------------------------------
//private member-functions

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
void SHA1::PadMessage()
{
	/*
	 *  Check to see if the current message block is too small to hold
	 *  the initial padding bits and length.  If so, we will pad the
	 *  block, process it, and then continue padding into a second
	 *  block.
	 */
	if (context.Message_Block_Index > 55)
	{
		context.Message_Block[context.Message_Block_Index++] = 0x80;
		while(context.Message_Block_Index < 64)
		{
			context.Message_Block[context.Message_Block_Index++] = 0;
		}

		ProcessMessageBlock();

		while(context.Message_Block_Index < 56)
		{
			context.Message_Block[context.Message_Block_Index++] = 0;
		}
	}
	else
	{
		context.Message_Block[context.Message_Block_Index++] = 0x80;
		while(context.Message_Block_Index < 56)
		{
			context.Message_Block[context.Message_Block_Index++] = 0;
		}
	}

	/*
	 *  Store the message length as the last 8 octets
	 */
	context.Message_Block[56] = context.Length_High >> 24;
	context.Message_Block[57] = context.Length_High >> 16;
	context.Message_Block[58] = context.Length_High >> 8;
	context.Message_Block[59] = context.Length_High;
	context.Message_Block[60] = context.Length_Low >> 24;
	context.Message_Block[61] = context.Length_Low >> 16;
	context.Message_Block[62] = context.Length_Low >> 8;
	context.Message_Block[63] = context.Length_Low;

	ProcessMessageBlock();
}

/**
 *  @brief      This member-function will process the next 512 bits of the
 *  		message stored in the Message_Block array.
 *
 *      	Many of the variable names in this code, especially the
 *      	single character names, were used because those were the
 *      	names used in the publication.
 *
 */  
void SHA1::ProcessMessageBlock()
{
	const uint32_t K[] =    {       /* Constants defined in SHA-1   */
		0x5A827999,
		0x6ED9EBA1,
		0x8F1BBCDC,
		0xCA62C1D6
	};
	int           t;                 /* Loop counter                */
	uint32_t      temp;              /* Temporary word value        */
	uint32_t      W[80];             /* Word sequence               */
	uint32_t      A, B, C, D, E;     /* Word buffers                */

	/*
	 *  Initialize the first 16 words in the array W
	 */
	for(t = 0; t < 16; t++)
	{
		W[t] = context.Message_Block[t * 4] << 24;
		W[t] |= context.Message_Block[t * 4 + 1] << 16;
		W[t] |= context.Message_Block[t * 4 + 2] << 8;
		W[t] |= context.Message_Block[t * 4 + 3];
	}

	for(t = 16; t < 80; t++)
	{
		W[t] = SHA1CircularShift(1,W[t-3] ^ W[t-8] ^ W[t-14] ^ W[t-16]);
	}

	A = context.Intermediate_Hash[0];
	B = context.Intermediate_Hash[1];
	C = context.Intermediate_Hash[2];
	D = context.Intermediate_Hash[3];
	E = context.Intermediate_Hash[4];

	for(t = 0; t < 20; t++)
	{
		temp =  SHA1CircularShift(5,A) +
			((B & C) | ((~B) & D)) + E + W[t] + K[0];
		E = D;
		D = C;
		C = SHA1CircularShift(30,B);
		B = A;
		A = temp;
	}

	for(t = 20; t < 40; t++)
	{
		temp = SHA1CircularShift(5,A) + (B ^ C ^ D) + E + W[t] + K[1];
		E = D;
		D = C;
		C = SHA1CircularShift(30,B);
		B = A;
		A = temp;
	}

	for(t = 40; t < 60; t++)
	{
		temp = SHA1CircularShift(5,A) +
			((B & C) | (B & D) | (C & D)) + E + W[t] + K[2];
		E = D;
		D = C;
		C = SHA1CircularShift(30,B);
		B = A;
		A = temp;
	}

	for(t = 60; t < 80; t++)
	{
		temp = SHA1CircularShift(5,A) + (B ^ C ^ D) + E + W[t] + K[3];
		E = D;
		D = C;
		C = SHA1CircularShift(30,B);
		B = A;
		A = temp;
	}

	context.Intermediate_Hash[0] += A;
	context.Intermediate_Hash[1] += B;
	context.Intermediate_Hash[2] += C;
	context.Intermediate_Hash[3] += D;
	context.Intermediate_Hash[4] += E;

	context.Message_Block_Index = 0;
}

//----------------------------------------------------------------------
//public member-functions

/**
 *  @brief 	Resets the sha1 context and starts a new
 *  		hashprocess
 *  @return	0 on succes an error number otherwise
 */  
int SHA1::Reset()
{

	context.Length_Low             = 0;
	context.Length_High            = 0;
	context.Message_Block_Index    = 0;

	context.Intermediate_Hash[0]   = 0x67452301;
	context.Intermediate_Hash[1]   = 0xEFCDAB89;
	context.Intermediate_Hash[2]   = 0x98BADCFE;
	context.Intermediate_Hash[3]   = 0x10325476;
	context.Intermediate_Hash[4]   = 0xC3D2E1F0;

	context.Computed   = 0;
	context.Corrupted  = 0;

	return shaSuccess;
}

/**
 *  @brief 	Data input.
 *
 *  		This memberfunction add data to the context.
 *
 *  @param	message_array The data to add
 *  @param	length The length of the data to add
 */  
int SHA1::Input(    const uint8_t  *message_array,
			size_t   length)
{
	if (!length)
	{
		return shaSuccess;
	}

	if (!message_array)
	{
		return shaNull;
	}

	if (context.Computed)
	{
		context.Corrupted = shaStateError;
		return shaStateError;
	}

	if (context.Corrupted)
	{
		return context.Corrupted;
	}
	while(length-- && !context.Corrupted)
	{
		context.Message_Block[context.Message_Block_Index++] =
			(*message_array & 0xFF);

		context.Length_Low += 8;
		if (context.Length_Low == 0)
		{
			context.Length_High++;
			if (context.Length_High == 0)
			{
				/* Message is too long */
				context.Corrupted = 1;
			}
		}

		if (context.Message_Block_Index == 64)
		{
			ProcessMessageBlock();
		}

		message_array++;
	}

	return shaSuccess;
}

/**
 *  @brief 	This ends the sha operation, zeroizing the context
 *  		and returning the computed hash.
 *
 *  @param	Message_Digest This is an OUT parameter which
 *  		contains the hash after the menberfunction returns
 *  @return	0 on succes, an error-code otherwise
 */  
int SHA1::Result( uint8_t 	  Message_Digest[HashSize])
{
	size_t i;

	if (!Message_Digest)
	{
		return shaNull;
	}

	if (context.Corrupted)
	{
		return context.Corrupted;
	}

	if (!context.Computed)
	{
		PadMessage();
		for(i=0; i<64; ++i)
		{
			/* message may be sensitive, clear it out */
			context.Message_Block[i] = 0;
		}
		context.Length_Low = 0;    /* and clear length */
		context.Length_High = 0;
		context.Computed = 1;
	}

	for(i = 0; i < HashSize; ++i)
	{
		Message_Digest[i] = context.Intermediate_Hash[i>>2]
			>> 8 * ( 3 - ( i & 0x03 ) );
	}

	return shaSuccess;
}


/**
 *  @brief 	This computes a sha1 digest on a string in one go.
 *
 *  @param	input_str the std::string input string
 *  @return	an std::string of HashSize bytes (20)
 */  
const std::string &fromString( const std::string &input_str, std::string &output_str) {
	SHA1 digest;
	uint8_t Message_Digest[HashSize];

	digest.Input((unsigned char*) input_str.c_str(), input_str.length());
	digest.Result(Message_Digest);
	output_str.assign ((char *)&Message_Digest[0],HashSize);
	return output_str;
}

/**
 *  @brief 	This computes a sha1 digest on a set of bytes in one go.
 *
 *  @param	bytes pointer to bytes
 *  @param	length number of bytes pointed to by the bytes pointer
 *  @return	an std::string of HashSize bytes (20)
 */  
const std::string &fromBytes( const unsigned char *bytes, const size_t length, std::string &output_str) {
	SHA1 digest;
	uint8_t Message_Digest[HashSize];

	digest.Input(bytes, length);
	digest.Result(Message_Digest);
	output_str.assign ((char *)&Message_Digest[0],HashSize);
	return output_str;
}


//----------------------------------------------------------------------
// end of namespace SHA1
}

//----------------------------------------------------------------------
//EOF
