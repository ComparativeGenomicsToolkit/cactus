#ifndef CACTUS_META_EVENT_H_
#define CACTUS_META_EVENT_H_

#include "cactusGlobals.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Meta event functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs a meta event, which contains all the essential info for a event.
 */
MetaEvent *metaEvent_construct(const char *header,
		CactusDisk *cactusDisk);

/*
 * Gets the name of the event.
 */
Name metaEvent_getName(MetaEvent *metaEvent);

/*
 * Gets the header line associated with the meta event.
 */
const char *metaEvent_getHeader(MetaEvent *metaEvent);

#endif
