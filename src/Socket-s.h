/*

   Socket.h		Socket functions

   This module contains an simple socket functions.

   Copyright (C) 2006, Wong Chi Kwong.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

#ifndef __SOCKET_H__
#define __SOCKET_H__

#define INET_SOCKET		0
#define LOCAL_SOCKET	1

#define SOCKET_MAX_PATH_LENGTH	64

#define SOCKET_REDIRECT_MAX_LENGTH	255

#define SOCKET_REDIRECT_STDOUT	'o'
#define SOCKET_REDIRECT_STDERR	'e'
#define SOCKET_REDIRECT_UNKNOWN	'u'

typedef struct Socket {
	int socketType;
	char *socketPath;
	int masterSocket;
	int currentSocket;
} Socket;

Socket *SocketCreate(const int socketType, const char *socketPath);
int SocketAcceptConnection(Socket* socketPointer);
Socket *SocketInitiateConnection(const int socketType, const char *socketPath);
int SocketReceive(Socket* socketPointer, void* address, const int maxNumOfByte);
int SocketSend(Socket* socketPointer, const void* address, const int numOfByte);
int SocketEndConnection(Socket* socketPointer);
int SocketFree(Socket* socketPointer);

void SocketSetRedirect(Socket* socketPointer);
void SocketResetRedirect();
int Socketfprintf(FILE *stream, const char * fmt,...);
int Socketprintf(const char * fmt,...);
int SocketRedirect();

#endif
