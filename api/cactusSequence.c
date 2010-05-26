#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic sequence functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

Sequence *sequence_construct(MetaSequence *metaSequence, Net *net) {
	Sequence *sequence;
	sequence = st_malloc(sizeof(Sequence));
	sequence->metaSequence = metaSequence;
	sequence->net = net;
	net_addSequence(net, sequence);
	return sequence;
}

void sequence_destruct(Sequence *sequence) {
	net_removeSequence(sequence_getNet(sequence), sequence);
	free(sequence);
}

MetaSequence *sequence_getMetaSequence(Sequence *sequence) {
	return sequence->metaSequence;
}

int32_t sequence_getStart(Sequence *sequence) {
	return metaSequence_getStart(sequence->metaSequence);
}

int32_t sequence_getLength(Sequence *sequence) {
	return metaSequence_getLength(sequence->metaSequence);
}

Name sequence_getName(Sequence *sequence) {
	return metaSequence_getName(sequence->metaSequence);
}

Event *sequence_getEvent(Sequence *sequence) {
	return eventTree_getEvent(net_getEventTree(sequence_getNet(sequence)), metaSequence_getEventName(sequence->metaSequence));
}

Net *sequence_getNet(Sequence *sequence) {
	return sequence->net;
}

char *sequence_getString(Sequence *sequence, int32_t start, int32_t length, int32_t strand) {
	return metaSequence_getString(sequence->metaSequence, start, length, strand);
}

const char *sequence_getHeader(Sequence *sequence) {
	return metaSequence_getHeader(sequence->metaSequence);
}

void sequence_check(Sequence *sequence) {
	Net *net = sequence_getNet(sequence);
	assert(net_getSequence(net, sequence_getName(sequence)) == sequence); //properly connected to the net..

	Group *parentGroup = net_getParentGroup(net);
	if(parentGroup != NULL) {
		Net *parentNet = group_getNet(parentGroup);
		Sequence *parentSequence = net_getSequence(parentNet, sequence_getName(sequence));
		if(parentSequence != NULL) {
			assert(event_getName(sequence_getEvent(sequence)) == event_getName(sequence_getEvent(parentSequence)));
		}
	}
}

/*
 * Serialisation functions.
 */

void sequence_writeBinaryRepresentation(Sequence *sequence, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	binaryRepresentation_writeElementType(CODE_SEQUENCE, writeFn);
	binaryRepresentation_writeName(sequence_getName(sequence), writeFn);
}

Sequence *sequence_loadFromBinaryRepresentation(void **binaryString, Net *net) {
	Sequence *sequence;

	sequence = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_SEQUENCE) {
		binaryRepresentation_popNextElementType(binaryString);
		sequence = sequence_construct(netDisk_getMetaSequence(net_getNetDisk(net), binaryRepresentation_getName(binaryString)), net);
	}
	return sequence;
}

Sequence *sequence_getStaticNameWrapper(Name name) {
	static Sequence sequence;
	static MetaSequence metaSequence;
	sequence.metaSequence = &metaSequence;
	metaSequence.name = name;
	return &sequence;
}
