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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "MemManager.h"
#include "Socket.h"

char socketRedirectTemp[SOCKET_REDIRECT_MAX_LENGTH+2];
Socket* socketRedirectPointer = NULL;


#ifdef WIN32
#else
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/un.h>
#endif

#ifdef WIN32

// WIN32 portion is a snub

Socket *SocketCreate(const int socketType, const char *socketPath) {
	return NULL;
}

int SocketAcceptConnection(Socket* socketPointer) {
	return -1;
}

Socket *SocketInitiateConnection(const int socketType, const char *socketPath) {
	return NULL;
}

int SocketReceive(Socket* socketPointer, void* address, const int maxNumOfByte) {
	return -1;
}

int SocketSend(Socket* socketPointer, const void* address, const int numOfByte) {
	return -1;
}

int SocketEndConnection(Socket* socketPointer) {
	return -1;
}

int SocketFree(Socket* socketPointer) {
	return -1;
}

#else


Socket *SocketCreate(const int socketType, const char *socketPath) {

	Socket *socketPointer;
	struct sockaddr_un local;
	int masterSocket;
	int length;

	if (socketType != LOCAL_SOCKET) {
		fprintf(stderr, "Only local socket is supported!\n");
		exit(1);
	}
	if (socketPath == '\0' || socketPath[0] == '\0') {
		fprintf(stderr, "Socket path is not specified!\n");
		exit(1);
	}

	masterSocket = socket(AF_UNIX, SOCK_STREAM, 0);
	if (masterSocket == -1) {
		return NULL;
	}

	local.sun_family = AF_UNIX;
	strncpy(local.sun_path, socketPath, SOCKET_MAX_PATH_LENGTH);
	local.sun_path[SOCKET_MAX_PATH_LENGTH] = '\0';
	unlink(local.sun_path);
	length = strlen(local.sun_path) + sizeof(local.sun_family);
	if (bind(masterSocket, (struct sockaddr *)&local, length) == -1) {
		return NULL;
	}

	socketPointer = MMUnitAllocate(sizeof(Socket));
	socketPointer->masterSocket = masterSocket;
	socketPointer->currentSocket = -1;
	socketPointer->socketType = socketType;
	socketPointer->socketPath = MMUnitAllocate(strlen(socketPath) + 1);
	strcpy(socketPointer->socketPath, socketPath);

	return socketPointer;

}

int SocketAcceptConnection(Socket* socketPointer) {

	struct sockaddr_un remote;
	int t;
	int n;

	n = listen(socketPointer->masterSocket, 0);
	if (n != 0) {
		return n;
	}

	t = sizeof(remote);
	socketPointer->currentSocket = accept(socketPointer->masterSocket, (struct sockaddr *)&remote, &t);
	if (socketPointer->currentSocket == -1) {
		return -1;
	}

	return 0;

}

Socket *SocketInitiateConnection(const int socketType, const char *socketPath) {

	int length;
	struct sockaddr_un remote;
	int currentSocket;
	Socket *socketPointer;

	if (socketType != LOCAL_SOCKET) {
		fprintf(stderr, "Only local socket is supported!\n");
		exit(1);
	}
	if (socketPath == '\0' || socketPath[0] == '\0') {
		fprintf(stderr, "Socket path is not specified!\n");
		exit(1);
	}

	currentSocket = socket(AF_UNIX, SOCK_STREAM, 0);
	if (currentSocket == -1) {
		return NULL;
	}

	remote.sun_family = AF_UNIX;
	strncpy(remote.sun_path, socketPath, SOCKET_MAX_PATH_LENGTH);
	length = strlen(remote.sun_path) + sizeof(remote.sun_family);

	if (connect(currentSocket, (struct sockaddr *)&remote, length) == -1) {
		return NULL;
	}

	socketPointer = MMUnitAllocate(sizeof(Socket));
	socketPointer->masterSocket = -1;
	socketPointer->currentSocket = currentSocket;
	socketPointer->socketType = socketType;
	socketPointer->socketPath = MMUnitAllocate(strlen(socketPath) + 1);
	strcpy(socketPointer->socketPath, socketPath);

	return socketPointer;

}

int SocketReceive(Socket* socketPointer, void *address, const int maxNumOfByte) {

	return recv(socketPointer->currentSocket, address, maxNumOfByte, 0);

}

int SocketSend(Socket* socketPointer, const void* address, const int numOfByte) {

	return send(socketPointer->currentSocket, address, numOfByte, 0);

}

int SocketEndConnection(Socket* socketPointer) {

	int n = 0;
	if (socketPointer->currentSocket > 0) {
		n = close(socketPointer->currentSocket);
	}
	socketPointer->currentSocket = -1;

	return n;

}

int SocketFree(Socket* socketPointer) {

	int n = 0;
	if (socketPointer->masterSocket > 0) {
		n = close(socketPointer->masterSocket);
	}
	MMUnitFree(socketPointer->socketPath, strlen(socketPointer->socketPath + 1));
	MMUnitFree(socketPointer, sizeof(Socket));

	return n;

}

#endif

void SocketSetRedirect(Socket* socketPointer) {

	socketRedirectPointer = socketPointer;

}

void SocketResetRedirect() {

	socketRedirectPointer = NULL;

}

int Socketfprintf(FILE *stream, const char * fmt,...) {

	va_list ap;
	int n;
    
	va_start(ap,fmt);

	if (socketRedirectPointer == NULL) {
		// Function as fprintf
		return vfprintf(stream, fmt, ap);
	} else {
		// Send through socket
		n = vsnprintf(socketRedirectTemp + 2, SOCKET_REDIRECT_MAX_LENGTH-1, fmt, ap);
		if (n < 0) {
			for (n=0; n<SOCKET_REDIRECT_MAX_LENGTH; n++) {
				socketRedirectTemp[n+2] = ' ';
			}
			vsnprintf(socketRedirectTemp + 2, SOCKET_REDIRECT_MAX_LENGTH-1, fmt, ap);
			n = SOCKET_REDIRECT_MAX_LENGTH;
		}
		if (stream == stdout) {
			socketRedirectTemp[0] = SOCKET_REDIRECT_STDOUT;
		} else if (stream == stderr) {
			socketRedirectTemp[0] = SOCKET_REDIRECT_STDERR;
		} else {
			socketRedirectTemp[0] = SOCKET_REDIRECT_UNKNOWN;
		}
		socketRedirectTemp[1] = (char)n;
		n = SocketSend(socketRedirectPointer, socketRedirectTemp, n+2);
		if (n < 0) {
			fprintf(stderr, "Cannot redirect through socket! The original text:\n");
			vfprintf(stderr, fmt, ap);
			return -1;
		}
	}

	return n;

}

int Socketprintf(const char * fmt,...) {

	va_list ap;
	int n;
    
	va_start(ap,fmt);

	if (socketRedirectPointer == NULL) {
		// Function as printf
		return vprintf(fmt, ap);
	} else {
		// Send through socket
		n = vsnprintf(socketRedirectTemp + 2, SOCKET_REDIRECT_MAX_LENGTH-1, fmt, ap);
		if (n < 0) {
			for (n=0; n<SOCKET_REDIRECT_MAX_LENGTH; n++) {
				socketRedirectTemp[n+2] = ' ';
			}
			vsnprintf(socketRedirectTemp + 2, SOCKET_REDIRECT_MAX_LENGTH-1, fmt, ap);
			n = SOCKET_REDIRECT_MAX_LENGTH;
		}
		socketRedirectTemp[0] = SOCKET_REDIRECT_STDOUT;
		socketRedirectTemp[1] = (char)n;
		n = SocketSend(socketRedirectPointer, socketRedirectTemp, n+2);
		if (n < 0) {
			fprintf(stderr, "Cannot redirect through socket! The original text:\n");
			vfprintf(stderr, fmt, ap);
			return -1;
		}
	}

	return n;

}

int SocketRedirect() {

	int n;
	int i;

	n = SocketReceive(socketRedirectPointer, socketRedirectTemp, 2);
	if (n <= 0) {
		return n;
	}
	if (n < 2) {
		fprintf(stderr, "Cannot receive length of text!\n");
		return -1;
	}

	n = (unsigned char)socketRedirectTemp[1];
	if (n == 0) {
		return 0;
	}

	n = SocketReceive(socketRedirectPointer, socketRedirectTemp+2, n);
	if (n <= 0) {
		return n;
	}

	if (socketRedirectTemp[0] == SOCKET_REDIRECT_STDOUT) {
		for (i=0; i<n; i++) {
			fprintf(stdout, "%c", socketRedirectTemp[i+2]);
		}
	} else if (socketRedirectTemp[0] == SOCKET_REDIRECT_STDERR) {
		for (i=0; i<n; i++) {
			fprintf(stderr, "%c", socketRedirectTemp[i+2]);
		}
	} else {
		// reserved
		fprintf(stderr, "Unknown stream type! The original text:\n");
		for (i=0; i<n; i++) {
			fprintf(stderr, "%c", socketRedirectTemp[i+2]);
		}
	}
	return n;

}

