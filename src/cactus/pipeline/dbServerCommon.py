#!/usr/bin/env python3

"""
Functions to get network info
"""

import platform
import socket
from contextlib import closing

from toil.lib.bioio import logger
from cactus.shared.common import cactus_call

def getHostName():
    if platform.system() == 'Darwin':
        # macOS doesn't have a true Docker bridging mode, so each
        # container is NATed with its own local IP and can't bind
        # to/connect to this computer's external IP. We have to
        # hardcode the loopback IP to work around this.
        return '127.0.0.1'
    return getPublicIP()

# Borrowed from toil code.
def getPublicIP():
    """Get the IP that this machine uses to contact the internet.
    If behind a NAT, this will still be this computer's IP, and not the router's."""
    try:
        # Try to get the internet-facing IP by attempting a connection
        # to a non-existent server and reading what IP was used.
        with closing(socket.socket(socket.AF_INET, socket.SOCK_DGRAM)) as sock:
            # 203.0.113.0/24 is reserved as TEST-NET-3 by RFC 5737, so
            # there is guaranteed to be no one listening on the other
            # end (and we won't accidentally DOS anyone).
            sock.connect(('203.0.113.1', 1))
            ip = sock.getsockname()[0]
        return ip
    except:
        # Something went terribly wrong. Just give loopback rather
        # than killing everything, because this is often called just
        # to provide a default argument
        return '127.0.0.1'

def findOccupiedPorts():
    """Attempt to find all currently taken TCP ports.

    Returns a set of ints, representing taken ports."""
    netstatOutput = cactus_call(parameters=["netstat", "-tuplen"], check_output=True)
    ports = set()
    for line in netstatOutput.split("\n"):
        fields = line.split()
        if len(fields) != 9:
            # Header or other garbage line
            continue
        port = int(fields[3].split(':')[-1])
        ports.add(port)
    logger.debug('Detected ports in use: %s' % repr(ports))
    return ports