#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>

#include "hashTableC.h"
#include "net2.h"
#include "net2Private.h"
#include "commonC.h"
#include "avl.h"

/*
 * Implementation of the net API
 *
 */

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic sequence functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

Sequence *sequence_construct(const char *name, int32_t length,
		const char *file, NetDisk *netDisk) {
	Sequence *sequence;
#ifdef BEN_DEBUG
	assert(length >= 0);
	assert(netDisk != NULL);
#endif
	sequence = malloc(sizeof(Sequence));
	sequence->file = stringCopy(file);
	sequence->name = stringCopy(name);
	sequence->length = length;
	sequence->netDisk = netDisk;

	netDisk_addSequence(netDisk, sequence);

	return sequence;
}

void sequence_destruct(Sequence *sequence) {
	netDisk_unloadSequence(sequence->netDisk, sequence);
	free(sequence->file);
	free(sequence->name);
	free(sequence);
}

int32_t sequence_getLength(Sequence *sequence) {
	return sequence->length;
}

const char *sequence_getName(Sequence *sequence) {
	return sequence->name;
}

const char *sequence_getFile(Sequence *sequence) {
	return sequence->file;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic end instance functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

EndInstance *endInstance_construct(const char *instance, End *end) {
	EndInstance *endInstance;
#ifdef BEN_DEBUG
	assert(insinstance != NULL);
#endif

	endInstance = malloc(sizeof(EndInstance));
	endInstance->instance = stringCopy(instance);
	endInstance->end = end;

	endInstance->coordinate = INT32_MAX;
	endInstance->sequence = NULL;
	endInstance->adjacency = NULL;
	endInstance->adjacency2 = NULL;
	endInstance->operation = NULL;
	endInstance->atomInstance = NULL;
	endInstance->parent = NULL;
	endInstance->children = constructEmptyList(0, NULL);

	end_addInstance(end, endInstance);
	return endInstance;
}

EndInstance *endInstance_constructWithCoordinates(const char *instance, End *end, int32_t coordinate, Sequence *sequence) {
	EndInstance *endInstance;
	endInstance = endInstance_construct(instance, end);
	endInstance->coordinate = coordinate;
	endInstance->sequence = sequence;
	return endInstance;
}

void endInstance_destruct(EndInstance *endInstance) {
	end_removeInstance(endInstance_getEnd(endInstance), endInstance);
	endInstance_breakAdjacency1(endInstance);
	endInstance_breakAdjacency2(endInstance);
	destructList(endInstance->children);
	free(endInstance->instance);
	free(endInstance);
}

const char *endInstance_getName(EndInstance *endInstance) {
	return endInstance->instance;
}

End *endInstance_getEnd(EndInstance *endInstance) {
	return endInstance->end;
}

AtomInstance *endInstance_getAtomInstance(EndInstance *endInstance) {
	return endInstance->atomInstance;
}

int32_t endInstance_getCoordinate(EndInstance *endInstance) {
	return endInstance->coordinate;
}

Sequence *endInstance_getSequence(EndInstance *endInstance) {
	return endInstance->sequence;
}

void endInstance_makeAdjacent1(EndInstance *endInstance, EndInstance *endInstance2) {
	endInstance_breakAdjacency1(endInstance);
	endInstance_breakAdjacency1(endInstance2);
	endInstance->adjacency = endInstance2;
	endInstance2->adjacency = endInstance;
}

void endInstance_makeAdjacent2(EndInstance *endInstance, EndInstance *endInstance2) {
	endInstance_breakAdjacency2(endInstance);
	endInstance_breakAdjacency2(endInstance2);
	endInstance->adjacency2 = endInstance2;
	endInstance2->adjacency2 = endInstance;
}

EndInstance *endInstance_getAdjacency(EndInstance *endInstance) {
	return endInstance->adjacency;
}

EndInstance *endInstance_getAdjacency2(EndInstance *endInstance) {
	return endInstance->adjacency2;
}

Operation *endInstance_getOperation(EndInstance *endInstance) {
	return endInstance->operation;
}

EndInstance *endInstance_getParent(EndInstance *endInstance) {
	return endInstance->parent;
}

int32_t endInstance_getChildNumber(EndInstance *endInstance) {
	return endInstance->children->length;
}

EndInstance *endInstance_getChild(EndInstance *endInstance, int32_t index) {
#ifdef BEN_DEBUG
	assert(endInstance_getChildNumber(endInstance) > index);
	assert(index >= 0);
#endif
	return endInstance->children->list[index];
}

void endInstance_makeParentAndChild(EndInstance *endInstanceParent, EndInstance *endInstanceChild) {
	listAppend(endInstanceParent->children, endInstanceChild);
	endInstanceChild->parent = endInstanceParent;
}

int32_t endInstance_isInternal(EndInstance *endInstance) {
	return endInstance_getChildNumber(endInstance) >= 0;
}

int32_t endInstance_isAugmented(EndInstance *endInstance) {

}

/*
 * Private functions.
 */

void endInstance_setAtomInstance(EndInstance *endInstance, AtomInstance *atomInstance) {
	endInstance->atomInstance = atomInstance; //done reciprocally by calling function.
}

void endInstance_setOperation(EndInstance *endInstance, Operation *operation) {
	endInstance->operation = operation;
}

void endInstance_breakAdjacency1(EndInstance *endInstance) {
	EndInstance *endInstance2;
	endInstance2 = endInstance_getAdjacency(endInstance);
	if(endInstance2 != NULL) {
		endInstance2->adjacency = NULL;
		endInstance->adjacency = NULL;
	}
}

void endInstance_breakAdjacency2(EndInstance *endInstance) {
	EndInstance *endInstance2;
	endInstance2 = endInstance_getAdjacency2(endInstance);
	if(endInstance2 != NULL) {
		endInstance2->adjacency2 = NULL;
		endInstance->adjacency2 = NULL;
	}
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic end functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int32_t end_constructP(const void *o1, const void *o2, void *a) {
	return strcmp(endInstance_getName((EndInstance *)o1), endInstance_getName((EndInstance *)o2));
}

End *end_construct(const char *name, Net *net) {
	End *end;
	end = malloc(sizeof(End));
	end->endInstances = avl_create(end_constructP, NULL, NULL);
	end->attachedAtom = NULL;
	net_addEnd(net, end);
	return end;
}

End *end_copyConstruct(End *end, Net *newNet) {
	End *end2 = end_construct(end_getName(end), newNet);
	End_InstanceIterator *iterator;
	EndInstance *endInstance;
	EndInstance *endInstance2;

	//Copy the instances.
	iterator = end_getInstanceIterator(end);
	while((endInstance = end_getNext(iterator)) != NULL) {
		if(endInstance_getCoordinate(endInstance) != INT32_MAX) {
			endInstance_constructWithCoordinates(endInstance_getName(endInstance), end2,
					endInstance_getCoordinate(endInstance), endInstance_getSequence(endInstance));
		}
		else {
			endInstance_construct(endInstance_getName(endInstance), end2);
		}
	}
	end_destructInstanceIterator(iterator);

	//Copy any parent child links.
	iterator = end_getInstanceIterator(end);
	while((endInstance = end_getNext(iterator)) != NULL) {
		if((endInstance2 = endInstance_getParent(endInstance)) != NULL) {
			endInstance_linkParentAndChild(end_getInstance(end2, endInstance_getName(endInstance2)),
										   end_getInstance(end2, endInstance_getName(endInstance)));
		}
	}
	end_destructInstanceIterator(iterator);
	return end2;
}

void end_destruct(End *end) {
	EndInstance *endInstance;
	net_removeEnd(end_getNet(end), end);
	while((endInstance = end_getFirst(end)) != NULL) {
		endInstance_destruct(endInstance);
	}
}

const char *end_getName(End *end) {
	return end->name;
}

Net *end_getNet(End *end) {
	return end->net;
}

AdjacencyComponent *end_getAdjacencyComponent(End *end) {
	return end->adjacencyComponent;
}

int32_t end_getInstanceNumber(End *end) {
	return avl_count(end->endInstances);
}

EndInstance *end_getInstance(End *end, const char *name) {
	static EndInstance endInstance;
	endInstance.instance = (char *)name;
	return avl_find(end->endInstances, &endInstance);
}

EndInstance *end_getFirst(End *end) {
	static struct avl_traverser iterator;
	avl_t_init(&iterator, end->endInstances);
	return avl_t_first(&iterator, end->endInstances);
}

End_InstanceIterator *end_getInstanceIterator(End *end) {
	End_InstanceIterator *iterator;
	iterator = mallocLocal(sizeof(End_InstanceIterator));
	avl_t_init(iterator, end->endInstances);
	return iterator;
}

EndInstance *end_getNext(End_InstanceIterator *iterator) {
	return avl_t_next(iterator);
}

void end_destructInstanceIterator(End_InstanceIterator *iterator) {
	free(iterator);
}

int32_t end_isStub(End *end) {
	return end_getName(end)[0] == '$';
}

int32_t end_isCap(End *end) {
	return !end_isStub(end) && end_getAtom(end) != NULL;
}

int32_t end_isAtomEnd(End *end) {
	return end_getAtom(end) != NULL;
}

/*
 * Private functions
 */

void end_addInstance(End *end, EndInstance *endInstance) {
	avl_insert(end->endInstances, endInstance);
}

void end_removeInstance(End *end, EndInstance *endInstance) {
	avl_delete(end->endInstances, endInstance);
}

void end_setAdjacencyComponent(End *end, AdjacencyComponent *adjacencyComponent) {
	end->adjacencyComponent = adjacencyComponent;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic atom instance functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

AtomInstance *atomInstance_constructP(const char *instance, Atom *atom) {
	AtomInstance *atomInstance;
	atomInstance = malloc(sizeof(AtomInstance));
	atomInstance->rInstance = malloc(sizeof(AtomInstance));
	atomInstance->rInstance->rInstance = atomInstance;
	atomInstance->atom = atom;
	atomInstance->rInstance->atom = atom_getReverse(atom);
	return atomInstance;
}

AtomInstance *atomInstance_construct(const char *instance, Atom *atom) {
	AtomInstance *atomInstance;
	atomInstance = atomInstance_constructP(instance, atom);
	atomInstance->leftEndInstance =
			endInstance_construct(instance, atom_getLeft(atom));
	endInstance_setAtomInstance(atomInstance->leftEndInstance, atomInstance);

	atomInstance->rInstance->leftEndInstance =
				endInstance_construct(instance, atom_getRight(atom));
	endInstance_setAtomInstance(atomInstance->rInstance->leftEndInstance, atomInstance->rInstance);
	atom_addInstance(atom, atomInstance);
	return atomInstance;
}

AtomInstance *atomInstance_constructWithCoordinates(const char *instance, Atom *atom,
		int32_t startCoordinate, Sequence *sequence) {
	AtomInstance *atomInstance;

#ifdef BEN_DEBUG
	assert(startCoordinate >= 0);
	assert(startCoordinate + atom_getLength(atom) <= sequence_getLength(sequence));
#endif

	atomInstance = atomInstance_constructP(instance, atom);
	atomInstance->leftEndInstance =
			endInstance_constructWithCoordinates(instance, atom_getLeft(atom),
					startCoordinate, sequence);
	endInstance_setAtomInstance(atomInstance->leftEndInstance, atomInstance);

	atomInstance->rInstance->leftEndInstance =
				endInstance_constructWithCoordinates(instance, atom_getRight(atom),
						-(startCoordinate + atom_getLength(atom)), sequence);
	endInstance_setAtomInstance(atomInstance->rInstance->leftEndInstance, atomInstance->rInstance);
	atom_addInstance(atom, atomInstance);
	return atomInstance;
}

void atomInstance_destruct(AtomInstance *atomInstance) {
	atom_removeInstance(atomInstance_getAtom(atomInstance), atomInstance);
	atomInstance_getLeft(atomInstance)->atomInstance = NULL;
	atomInstance_getRight(atomInstance)->atomInstance = NULL;
	free(atomInstance->rInstance);
	free(atomInstance);
}

Atom *atomInstance_getAtom(AtomInstance *atomInstance) {
	return atomInstance->atom;
}

const char *atomInstance_getName(AtomInstance *atomInstance) {
	return endInstance_getName(atomInstance_getLeft(atomInstance));
}

AtomInstance *atomInstance_getReverse(AtomInstance *atomInstance) {
	return atomInstance->rInstance;
}

int32_t atomInstance_getStart(AtomInstance *atomInstance) {
	return endInstance_getCoordinate(atomInstance_getLeft(atomInstance));
}

int32_t atomInstance_getLength(AtomInstance *atomInstance) {
	return atom_getLength(atomInstance_getAtom(atomInstance));
}

Sequence *atomInstance_getSequence(AtomInstance *atomInstance) {
	return endInstance_getSequence(atomInstance_getLeft(atomInstance));
}

EndInstance *atomInstance_getLeft(AtomInstance *atomInstance) {
	return atomInstance->leftEndInstance;
}

EndInstance *atomInstance_getRight(AtomInstance *atomInstance) {
	return atomInstance_getLeft(atomInstance_getReverse(atomInstance));
}

AtomInstance *atomInstance_getParent(AtomInstance *atomInstance) {
	EndInstance *endInstance;
	AtomInstance *atomInstance2;
	endInstance = atomInstance_getLeft(atomInstance);
	while((endInstance = endInstance_getParent(endInstance)) != NULL) {
		if((atomInstance2 = endInstance_getAtomInstance(endInstance)) != NULL) {
			return atomInstance2;
		}
	}
	return NULL;
}

int32_t atomInstance_getChildNumber(AtomInstance *atomInstance) {
	return 0;
}

AtomInstance *atomInstance_getChild(AtomInstance *atomInstance, int32_t index) {
	return NULL;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic atom functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int32_t atomConstruct_constructP(const void *o1, const void *o2, void *a) {
	return strcmp(atomInstance_getName((AtomInstance *)o1), atomInstance_getName((AtomInstance *)o2));
}

Atom *atom_construct(const char *name, int32_t length, Net *net) {
	Atom *atom;
	atom = malloc(sizeof(Atom));
	atom->rAtom = malloc(sizeof(Atom));
	atom->rAtom->rAtom = atom;
	atom->atomContents = malloc(sizeof(struct AtomContents));
	atom->rAtom->atomContents = atom->atomContents;

	atom->atomContents->name = stringCopy(name);
	atom->atomContents->length = length;
	atom->atomContents->net = net;
	atom->atomContents->atomInstances = avl_create(atomConstruct_constructP, NULL, NULL);

	return atom;
}

void atom_destructP(AtomInstance *atomInstance, void *o) {
	atomInstance_destruct(atomInstance);
}

void atom_destruct(Atom *atom) {
	net_removeAtom(atom_getNet(atom), atom);
	free(atom->rAtom);
	free(atom->atomContents->name);
	avl_destroy(atom->atomContents->atomInstances, (void (*) (void *avl_item, void *avl_param))atom_destructP);
	free(atom->atomContents);
	free(atom);
}

const char *atom_getName(Atom *atom) {
	return atom->atomContents->name;
}

int32_t atom_getLength(Atom *atom) {
	return atom->atomContents->length;
}

Net *atom_getNet(Atom *atom) {
	return atom->atomContents->net;
}

End *atom_getLeft(Atom *atom) {
	return atom->leftEnd;
}

End *atom_getRight(Atom *atom) {
	return atom->rAtom->leftEnd;
}

Atom *atom_getReverse(Atom *atom) {
	return atom->rAtom;
}

int32_t atom_getInstanceNumber(Atom *atom) {
	return avl_count(atom->atomContents->atomInstances);
}

AtomInstance *atom_getInstance(Atom *atom, const char *name) {
	static AtomInstance atomInstance;
	//atomInstance.instance
	//return avl_find(atom->atomContents->atomInstances, );
}

AtomInstance *atom_getFirst(Atom *atom) {
	//return
}

Atom_InstanceIterator *atom_getInstanceIterator(Atom *atom) {
	Atom_InstanceIterator *iterator;
	iterator = malloc(sizeof(Atom_InstanceIterator));
	avl_t_init(iterator, atom->atomContents->atomInstances);
	return iterator;
}

AtomInstance *atom_getNext(Atom_InstanceIterator *iterator) {
	return avl_t_next(iterator);
}

void atom_destructInstanceIterator(Atom_InstanceIterator *atomInstanceIterator) {
	free(atomInstanceIterator);
}

/*
 * Private functions.
 */

void atom_addInstance(Atom *atom, AtomInstance *atomInstance) {
	avl_insert(atom->atomContents->atomInstances, atomInstance);
}

void atom_removeInstance(Atom *atom, AtomInstance *atomInstance) {
	avl_delete(atom->atomContents->atomInstances, atomInstance);
}


////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic adjacency component functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int32_t adjacencyComponent_constructP(const void *o1, const void *o2, void *a) {
	return strcmp(end_getName((End *)o1), end_getName((End *)o2));
}

AdjacencyComponent *adjacencyComponent_construct(Net *net, Net *nestedNet) {
	AdjacencyComponent *adjacencyComponent;
	adjacencyComponent = malloc(sizeof(AdjacencyComponent));

	adjacencyComponent->net = net;
	adjacencyComponent->nestedNetName = stringCopy(net_getName(nestedNet));
	adjacencyComponent->ends = avl_create(adjacencyComponent_constructP, NULL, NULL);

	net_addAdjacencyComponent(net, adjacencyComponent);

	return adjacencyComponent;
}

void adjacencyComponent_destruct(AdjacencyComponent *adjacencyComponent) {
	net_removeChain(adjacencyComponent_getNet(adjacencyComponent), adjacencyComponent);
	net_setParentAdjacencyComponent(adjacencyComponent_getNestedNet(adjacencyComponent), NULL);
	free(adjacencyComponent->nestedNetName);
	free(adjacencyComponent);
}

Net *adjacencyComponent_getNet(AdjacencyComponent *adjacencyComponent) {
	return adjacencyComponent->net;
}

const char *adjacencyComponent_getNestedNetName(AdjacencyComponent *adjacencyComponent) {
	return adjacencyComponent->nestedNetName;
}

Net *adjacencyComponent_getNestedNet(AdjacencyComponent *adjacencyComponent) {
	return netDisk_getNet(net_getNetDisk(adjacencyComponent_getNestedNet(adjacencyComponent)), adjacencyComponent->nestedNetName);
}

Chain *adjacencyComponent_getChain(AdjacencyComponent *adjacencyComponent) {
	return adjacencyComponent->chain;
}

End *adjacencyComponent_getEnd(AdjacencyComponent *adjacencyComponent, const char *name) {
	static End end;
	end.name = (char *)name;
	return avl_find(adjacencyComponent->ends, &end);
}

int32_t adjacencyComponent_getEndNumber(AdjacencyComponent *adjacencyComponent) {
	return avl_count(adjacencyComponent->ends);
}

AdjacencyComponent_EndIterator *adjacencyComponent_getEndIterator(AdjacencyComponent *adjacencyComponent) {
	AdjacencyComponent_EndIterator *iterator;
	iterator = malloc(sizeof(AdjacencyComponent_EndIterator));
	avl_t_init(iterator, adjacencyComponent->ends);
	return iterator;
}

End *adjacencyComponent_getNextEnd(AdjacencyComponent_EndIterator *endIterator) {
	return avl_t_next(endIterator);
}

void adjacencyComponent_destructEndIterator(AdjacencyComponent_EndIterator *endIterator) {
	free(endIterator);
}

/*
 * Private functions.
 */

void adjacencyComponent_setChain(AdjacencyComponent *adjacencyComponent, Chain *chain) {
	adjacencyComponent->chain = chain;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic chain functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

Chain *chain_construct(End *leftEnd, End *rightEnd, AdjacencyComponent *adjacencyComponent, Net *net, Chain *pLink) {
	Chain *chain;
	chain = malloc(sizeof(Chain));

	chain->nLink = NULL;
	chain->pLink = pLink;
	if(pLink != NULL) {
		pLink->nLink = chain;
	}
	chain->adjacencyComponent = adjacencyComponent;
	chain->leftEnd = leftEnd;
	chain->rightEnd = rightEnd;
	chain->net = net;

	net_addChain(net, chain); //will set the chain indices.
	return chain;
}

void chain_destruct(Chain *chain, int32_t recursive) {
	Chain *nChain;
	if(chain) {
		net_removeChain(chain_getNet(chain), chain);
	}
	while(recursive && (chain = chain_getNextLink(chain)) != NULL) {
		chain_destruct(nChain, recursive);
	}
	free(chain);
}

Chain *chain_getNextLink(Chain *chain) {
	return chain->nLink;
}

Chain *chain_getPreviousLink(Chain *chain) {
	return chain->pLink;
}

AdjacencyComponent *chain_getAdjacencyComponent(Chain *chain) {
	return chain->adjacencyComponent;
}

End *chain_getLeft(Chain *chain) {
	return chain->leftEnd;
}

End *chain_getRight(Chain *chain) {
	return chain->rightEnd;
}

Net *chain_getNet(Chain *chain) {
	return chain->net;
}

int32_t chain_getIndex(Chain *chain) {
	return chain->chainIndex;
}

/*
 * Private functions
 */

void chain_setIndex(Chain *chain, int32_t index) {
	chain->chainIndex = index;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic operation functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs an operation.
 */
Operation *operation_construct(Net *net) {
	Operation *operation;
	operation = malloc(sizeof(Operation));

	net_addOperation(net, operation);
	return operation;
}

/*
 * Destructs the operation.
 */
void operation_destruct(Operation *operation) {
	free(operation);
}

/*
 * Gets the net it is part of.
 */
Net *operation_getNet(Operation *operation) {
	return operation->net;
}

/*
 * Private functions
 */
void operation_setIndex(Operation *operation, int32_t index) {
	operation->index = index;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic net functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int32_t net_constructSequencesP(const void *o1, const void *o2, void *a) {
	return strcmp(sequence_getName((Sequence *)o1), sequence_getName((Sequence *)o2));
}

int32_t net_constructEndsP(const void *o1, const void *o2, void *a) {
	return strcmp(end_getName((End *)o1), end_getName((End *)o2));
}

int32_t net_constructAtomsP(const void *o1, const void *o2, void *a) {
	return strcmp(atom_getName((Atom *)o1), atom_getName((Atom *)o2));
}

int32_t net_constructAdjacencyComponentsP(const void *o1, const void *o2, void *a) {
	return strcmp(adjacencyComponent_getNestedNetName((AdjacencyComponent *)o1),
			adjacencyComponent_getNestedNetName((AdjacencyComponent *)o2));
}

int32_t net_constructChainsP(const void *o1, const void *o2, void *a) {
	return chain_getIndex((Chain *)o1) - chain_getIndex((Chain *)o2);
}

int32_t net_constructOperationsP(const void *o1, const void *o2, void *a) {
	return operation_getIndex((Operation *)o1) - operation_getIndex((Operation *)o2);
}

Net *net_construct(const char *name, NetDisk *netDisk) {
	Net *net;
	net = malloc(sizeof(Net));

	net->name = stringCopy(name);

	net->sequences = avl_create(net_constructSequencesP, NULL, NULL);
	net->ends = avl_create(net_constructEndsP, NULL, NULL);
	net->atoms = avl_create(net_constructAtomsP, NULL, NULL);
	net->adjacencyComponents = avl_create(net_constructAdjacencyComponentsP, NULL, NULL);
	net->chains = avl_create(net_constructChainsP, NULL, NULL);
	net->operations = avl_create(net_constructOperationsP, NULL, NULL);

	net->parentNetName = NULL;
	net->netDisk = netDisk;

	netDisk_addNet(net->netDisk, net);

	return net;
}

const char *net_getName(Net *net) {
	return net->name;
}

NetDisk *net_getNetDisk(Net *net) {
	return net->netDisk;
}

void net_destructEndP(End *end, void *o) {
	end_destruct(end);
}

void net_destructAtomP(Atom *atom, void *o) {
	atom_destruct(atom);
}

void net_destructAdjacencyComponentP(AdjacencyComponent *adjacencyComponent, void *o) {
	adjacencyComponent_destruct(adjacencyComponent);
}

void net_destructChainP(Chain *chain, void *o) {
	chain_destruct(chain, TRUE);
}

void net_destructOperationP(Operation *operation, void *o) {
	operation_destruct(operation);
}

void net_destruct(Net *net, int32_t recursive) {
	Net_AdjacencyComponentIterator *iterator;
	AdjacencyComponent *adjacencyComponent;

	if(recursive) {
		iterator = net_getAdjacencyComponentIterator(net);
		while((adjacencyComponent = net_getNextAdjacencyComponent(iterator)) != NULL) {
			net_destruct(adjacencyComponent_getNestedNet(adjacencyComponent), recursive);
		}
		net_destructAdjacencyComponentIterator(iterator);
	}

	netDisk_unloadNet(net->netDisk, net);
	avl_destroy(net->sequences, NULL);
	avl_destroy(net->ends, (void (*) (void *avl_item, void *avl_param))net_destructEndP);
	avl_destroy(net->atoms, (void (*) (void *avl_item, void *avl_param))net_destructAtomP);
	avl_destroy(net->adjacencyComponents, (void (*) (void *avl_item, void *avl_param))net_destructAdjacencyComponentP);
	avl_destroy(net->chains, (void (*) (void *avl_item, void *avl_param))net_destructChainP);
	avl_destroy(net->operations, (void (*) (void *avl_item, void *avl_param))net_destructOperationP);

	free(net->name);
	if(net->parentNetName != NULL) {
		free(net->parentNetName);
	}
	free(net);
}

void net_addSequence(Net *net, Sequence *sequence) {
	avl_insert(net->sequences, sequence);
}

Sequence *net_getSequence(Net *net, const char *name) {
	static Sequence sequence;
	sequence.name = (char *)name;
	return avl_find(net->sequences, &sequence);
}

int32_t net_getSequenceNumber(Net *net) {
	return avl_count(net->sequences);
}

Net_SequenceIterator *net_getSequenceIterator(Net *net) {
	Net_SequenceIterator *iterator;
	iterator = malloc(sizeof(Net_SequenceIterator));
	avl_t_init(iterator, net->sequences);
	return iterator;
}

Sequence *net_getNextSequence(Net_SequenceIterator *sequenceIterator) {
	return avl_t_next(sequenceIterator);
}

void net_destructSequenceIterator(Net_SequenceIterator *sequenceIterator) {
	free(sequenceIterator);
}

End *net_getEnd(Net *net, const char *name) {
	static End end;
	end.name = (char *)name;
	return avl_find(net->ends, &end);
}

int32_t net_getEndNumber(Net *net) {
	return avl_count(net->ends);
}

Net_EndIterator *net_getEndIterator(Net *net) {
	Net_EndIterator *iterator;
	iterator = malloc(sizeof(Net_EndIterator));
	avl_t_init(iterator, net->ends);
	return iterator;
}

End *net_getNextEnd(Net_EndIterator *endIterator) {
	return avl_t_next(endIterator);
}

void net_destructEndIterator(Net_EndIterator *endIterator) {
	free(endIterator);
}

Atom *net_getAtom(Net *net, const char *name) {
	static Atom atom;
	atom.name = (char *)name;
	return avl_find(net->atoms, &atom);
}

int32_t net_getAtomNumber(Net *net) {
	return avl_count(net->atoms);
}

Net_AtomIterator *net_getAtomIterator(Net *net) {
	Net_SequenceIterator *iterator;
	iterator = malloc(sizeof(Net_AtomIterator));
	avl_t_init(iterator, net->atoms);
	return iterator;
}

Atom *net_getNextAtom(Net_AtomIterator *atomIterator) {
	return avl_t_next(atomIterator);
}

void net_destructAtomIterator(Net_AtomIterator *atomIterator) {
	free(atomIterator);
}

AdjacencyComponent *net_getAdjacencyComponent(Net *net, const char *netName) {
	static AdjacencyComponent adjacencyComponent;
	adjacencyComponent.nestedNetName = (char *)netName;
	return avl_find(net->adjacencyComponents, &adjacencyComponent);
}

int32_t net_getAdjacencyComponentNumber(Net *net) {
	return avl_count(net->adjacencyComponents);
}

Net_AdjacencyComponentIterator *net_getAdjacencyComponentIterator(Net *net) {
	Net_AdjacencyComponentIterator *iterator;
	iterator = malloc(sizeof(Net_AdjacencyComponentIterator));
	avl_t_init(iterator, net->adjacencyComponents);
	return iterator;
}

AdjacencyComponent *net_getNextAdjacencyComponent(Net_AdjacencyComponentIterator *adjacencyComponentIterator) {
	return avl_t_next(adjacencyComponentIterator);
}

void net_destructAdjacencyComponentIterator(Net_AdjacencyComponentIterator *adjacencyComponentIterator) {
	free(adjacencyComponentIterator);
}

AdjacencyComponent *net_getParentAdjacencyComponent(Net *net) {
	Net *net2;
	net2 = netDisk_getNet(net_getNetDisk(net), net->parentNetName);
	return net_getAdjacencyComponent(net2, net_getName(net));
}

Chain *net_getChain(Net *net, int32_t index) {
	static Chain chain;
	chain.chainIndex = index;
	return avl_find(net->chains, &chain);
}

int32_t net_getChainNumber(Net *net) {
	return avl_count(net->chains);
}

Net_ChainIterator *net_getChainIterator(Net *net) {
	Net_ChainIterator *iterator;
	iterator = malloc(sizeof(Net_ChainIterator));
	avl_t_init(iterator, net->chains);
	return iterator;
}

Chain *net_getNextChain(Net_ChainIterator *chainIterator) {
	return avl_t_next(chainIterator);
}

void net_destructChainIterator(Net_ChainIterator *chainIterator) {
	free(chainIterator);
}

Operation *net_getOperation(Net *net, int32_t index) {
	static Operation operation;
	operation.index = index;
	return avl_find(net->operations, &operation);
}

int32_t net_getOperationNumber(Net *net) {
	return avl_count(net->operations);
}

Net_OperationIterator *net_getOperationIterator(Net *net) {
	Net_OperationIterator *iterator;
	iterator = malloc(sizeof(Net_OperationIterator));
	avl_t_init(iterator, net->operations);
	return iterator;
}

Operation *net_getNextOperation(Net_OperationIterator *operationIterator) {
	return avl_t_next(operationIterator);
}

void net_destructOperationIterator(Net_OperationIterator *operationIterator) {
	free(operationIterator);
}

char *net_makeBinaryRepresentation(Net *net) {
}

char *net_makeXMLRepresentation(Net *net) {
}

/*
 * Private functions
 */

void net_addAtom(Net *net, Atom *atom) {
	avl_insert(net->atoms, atom);
}

void net_removeAtom(Net *net, Atom *atom) {
	avl_delete(net->atoms, atom);
}

void net_addEnd(Net *net, End *end) {
	avl_insert(net->ends, end);
}

void net_removeEnd(Net *net, End *end) {
	avl_delete(net->ends, end);
}

void net_addChain(Net *net, Chain *chain) {
	chain_setIndex(chain, net->chainIndex++);
	avl_insert(net->chains, chain);
}

void net_removeChain(Net *net, Chain *chain) {
	avl_delete(net->chains, chain);
}

void net_addAdjacencyComponent(Net *net, AdjacencyComponent *adjacencyComponent) {
	avl_insert(net->adjacencyComponents, adjacencyComponent);
}

void net_removeAdjacencyComponent(Net *net, AdjacencyComponent *adjacencyComponent) {
	avl_delete(net->adjacencyComponents, adjacencyComponent);
}

void net_setParentAdjacencyComponent(Net *net, AdjacencyComponent *adjacencyComponent) {
	net->parentNetName = stringCopy(net_getName(adjacencyComponent_getNet(adjacencyComponent)));
}

void net_addOperation(Net *net, Operation *operation) {
	operation_setIndex(operation, net_getOperationNumber(net));
	avl_insert(net->operations, operation);
}

void net_removeOperation(Net *net, Operation *operation) {
	avl_delete(net->operations, operation);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic net disk functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

NetDisk *netDisk_construct(const char *netDiskFile) {
}

void netDisk_destruct(NetDisk *netDisk) {
}

int32_t netDisk_write(NetDisk *netDisk) {

}

Net *netDisk_getNet(NetDisk *netDisk, const char *netName) {
}

Sequence *netDisk_getSequence(NetDisk *netDisk, const char *sequenceName) {
}

/*
 * Private functions.
 */

void netDisk_addSequence(NetDisk *netDisk, Net *net);

void netDisk_unloadSequence(NetDisk *netDisk, const char *sequenceName);

int32_t netDisk_netIsLoaded(NetDisk *netDisk, Net *net);

void netDisk_addNet(NetDisk *netDisk, Net *net);

void netDisk_unloadNet(NetDisk *netDisk, Net *net);

