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
//Misc reused functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void *getFirst(struct avl_table *items) {
	static struct avl_traverser iterator;
	avl_t_init(&iterator, items);
	return avl_t_first(&iterator, items);
}

struct avl_traverser *getIterator(struct avl_table *items) {
	struct avl_traverser *iterator;
	iterator = mallocLocal(sizeof(struct avl_traverser));
	avl_t_init(iterator, items);
	return iterator;
}

void destructIterator(struct avl_traverser *iterator) {
	free(iterator);
}

void *getNext(struct avl_traverser *iterator) {
	return avl_t_next(iterator);
}

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
	int32_t i;
	EndInstance *endInstance2;

	//Remove from end.
	end_removeInstance(endInstance_getEnd(endInstance), endInstance);
	//Break adjacencies.
	endInstance_breakAdjacency1(endInstance);
	endInstance_breakAdjacency2(endInstance);
	//Deal with parental and child references to this instance.
	if((endInstance2 = endInstance_getParent(endInstance)) != NULL) {
		listRemove(endInstance2->children, endInstance);
	}
	for(i=0; i<endInstance->children->length; i++) {
		endInstance2 = endInstance->children->list[i];
		endInstance2->parent = NULL;
	}
	destructList(endInstance->children);
	//Free name string
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
	if(!listContains(endInstanceParent->children, endInstanceChild)) { //defensive, means second calls will have no effect.
		listAppend(endInstanceParent->children, endInstanceChild);
	}
	endInstanceChild->parent = endInstanceParent;
}

int32_t endInstance_isInternal(EndInstance *endInstance) {
	return endInstance_getChildNumber(endInstance) >= 0;
}

int32_t endInstance_isAugmented(EndInstance *endInstance) {
	return end_getAtom(endInstance_getEnd(endInstance)) != NULL && endInstance_getAtomInstance(endInstance) == NULL;
}

/*
 * Private functions.
 */

void endInstance_setAtomInstance(EndInstance *endInstance, AtomInstance *atomInstance) {
	endInstance->atomInstance = atomInstance;
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
	end->adjacencyComponent = NULL;
	end->net = net;
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
	//remove from net.
	net_removeEnd(end_getNet(end), end);

	//remove instances
	while((endInstance = end_getFirst(end)) != NULL) {
		endInstance_destruct(endInstance);
	}
	//now the actual instances.
	avl_destroy(end->endInstances, NULL);

	free(end->name);
	free(end);
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
	return getFirst(end->endInstances);
}

End_InstanceIterator *end_getInstanceIterator(End *end) {
	return getIterator(end->endInstances);
}

EndInstance *end_getNext(End_InstanceIterator *iterator) {
	return getNext(iterator);
}

void end_destructInstanceIterator(End_InstanceIterator *iterator) {
	destructIterator(iterator);
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
	//argument may be NULL
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
	//remove from ends.
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
	atom->atomContents->atomInstances = avl_create(atomConstruct_constructP, NULL, NULL);
	atom->atomContents->length = length;
	atom->atomContents->net = net;

	return atom;
}

void atom_destruct(Atom *atom) {
	AtomInstance *atomInstance;
	//remove from net.
	net_removeAtom(atom_getNet(atom), atom);

	//remove instances
	while((atomInstance = atom_getFirst(atom)) != NULL) {
		atomInstance_destruct(atomInstance);
	}
	//now the actual instances.
	avl_destroy(atom->atomContents->atomInstances, NULL);

	free(atom->rAtom);
	free(atom->atomContents->name);
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
	static EndInstance endInstance;
	endInstance.instance = (char *)name;
	atomInstance.leftEndInstance = &endInstance;
	return avl_find(atom->atomContents->atomInstances, &atomInstance);
}

AtomInstance *atom_getFirst(Atom *atom) {
	return getFirst(atom->atomContents->atomInstances);
}

Atom_InstanceIterator *atom_getInstanceIterator(Atom *atom) {
	return getIterator(atom->atomContents->atomInstances);
}

AtomInstance *atom_getNext(Atom_InstanceIterator *iterator) {
	return getNext(iterator);
}

void atom_destructInstanceIterator(Atom_InstanceIterator *atomInstanceIterator) {
	destructIterator(atomInstanceIterator);
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
	adjacencyComponent->chain = NULL;
	adjacencyComponent->nestedNetName = stringCopy(net_getName(nestedNet));
	adjacencyComponent->ends = avl_create(adjacencyComponent_constructP, NULL, NULL);
	adjacencyComponent_updateContainedEnds(adjacencyComponent);

	net_addAdjacencyComponent(net, adjacencyComponent);

	return adjacencyComponent;
}

void adjacencyComponent_destructP(End *end, void *o) {
	end_setAdjacencyComponent(end, NULL);
}

void adjacencyComponent_updateContainedEnds(AdjacencyComponent *adjacencyComponent) {
	Net *net;
	Net_EndIterator *iterator;
	End *end;
	End *end2;
	//wipe the slate clean.
	avl_destroy(adjacencyComponent->ends, (void (*)(void *, void *))adjacencyComponent_destructP);
	adjacencyComponent->ends = avl_create(adjacencyComponent_constructP, NULL, NULL);
	//now calculate the ends
	net = adjacencyComponent_getNet(adjacencyComponent);
	iterator = net_getEndIterator(adjacencyComponent_getNestedNet(adjacencyComponent));
	while((end = net_getNextEnd(iterator)) != NULL) {
		if((end2 = net_getEnd(net, end_getName(end))) != NULL) {
			avl_insert(adjacencyComponent->ends, end2);
			end_setAdjacencyComponent(end2, adjacencyComponent);
		}
	}
	net_destructEndIterator(iterator);
}

void adjacencyComponent_destruct(AdjacencyComponent *adjacencyComponent) {
	//Detatch from the nets.
	net_removeAdjacencyComponent(adjacencyComponent_getNet(adjacencyComponent), adjacencyComponent);
	net_setParentAdjacencyComponent(adjacencyComponent_getNestedNet(adjacencyComponent), NULL);
	//Detatch the ends from their adjacency components.
	avl_destroy(adjacencyComponent->ends, (void (*)(void *, void *))adjacencyComponent_destructP);
	//Free the memory
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
	return getIterator(adjacencyComponent->ends);
}

End *adjacencyComponent_getNextEnd(AdjacencyComponent_EndIterator *endIterator) {
	return getNext(endIterator);
}

void adjacencyComponent_destructEndIterator(AdjacencyComponent_EndIterator *endIterator) {
	destructIterator(endIterator);
}

/*
 * Private functions.
 */

void adjacencyComponent_setChain(AdjacencyComponent *adjacencyComponent, Chain *chain) {
	//argument may be NULL
	adjacencyComponent->chain = chain;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic link functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

Link *link_construct(End *leftEnd, End *rightEnd, AdjacencyComponent *adjacencyComponent, Chain *parentChain) {
	Link *link;
	link = malloc(sizeof(Link));

	link->leftEnd = leftEnd;
	link->rightEnd = rightEnd;
	link->chain = parentChain;
	link->adjacencyComponent = adjacencyComponent;

	chain_addLink(parentChain, link); //will set the link indices.
	return link;
}

void link_destruct(Link *link) {
	Link *link2;
	link_getChain(link)->linkNumber = link->linkIndex;
	if(link->pLink == NULL) {
		link_getChain(link)->link = NULL;
	}
	else {
		link->pLink->nLink = NULL;
	}
	while(link != NULL) {
		link2 = link->nLink;
		link = link->nLink;
		free(link2);
	}
}

Link *link_getNextLink(Link *link) {
	return link->nLink;
}

Link *link_getPreviousLink(Link *link) {
	return link->pLink;
}

AdjacencyComponent *link_getAdjacencyComponent(Link *link) {
	return link->adjacencyComponent;
}

End *link_getLeft(Link *link) {
	return link->leftEnd;
}

End *link_getRight(Link *link) {
	return link->rightEnd;
}

Chain *link_getChain(Link *link) {
	return link->chain;
}

int32_t link_getIndex(Link *link) {
	return link->linkIndex;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic chain functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

Chain *chain_construct(Net *net) {
	Chain *chain;
	chain = malloc(sizeof(Chain));
	chain->net = net;
	chain->link = NULL;
	chain->linkNumber = 0;
	net_addChain(net, chain);
	return chain;
}

void chain_destruct(Chain *chain) {
	net_removeChain(chain_getNet(chain), chain);
	if(chain->link != NULL) {
		link_destruct(chain->link);
	}
	free(chain);
}

Link *chain_getLink(Chain *chain, int32_t linkIndex) {
	int32_t i;
	Link *link;

#ifdef BEN_DEBUG
	assert(linkIndex >= 0);
	assert(linkIndex < chain->linkNumber);
#endif

	i=0;
	link = chain->link;
	while(i++ < linkIndex) {
		link = link->nLink;
	}
	return link;
}

int32_t chain_getLength(Chain *chain) {
	return chain->linkNumber;
}

int32_t chain_getIndex(Chain *chain) {
	return chain->chainIndex;
}

Net *chain_getNet(Chain *chain) {
	return chain->net;
}

/*
 * Private functions
 */

void chain_addLink(Chain *chain, Link *childLink) {
	Link *pLink;
	if(chain->linkNumber != 0) {
		pLink = chain_getLink(chain, chain->linkNumber -1);
		pLink->nLink = childLink;
		childLink->pLink = pLink;
	}
	else {
		childLink->pLink = NULL;
		chain->link = childLink;
	}
	childLink->nLink = NULL;
	childLink->linkIndex = chain->linkNumber++;
}

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
	net->operationIndex = 0;
	net->chainIndex = 0;

	netDisk_addNet(net->netDisk, net);

	return net;
}

const char *net_getName(Net *net) {
	return net->name;
}

NetDisk *net_getNetDisk(Net *net) {
	return net->netDisk;
}

void net_destruct(Net *net, int32_t recursive) {
	Net_AdjacencyComponentIterator *iterator;
	End *end;
	Atom *atom;
	AdjacencyComponent *adjacencyComponent;
	Chain *chain;
	Operation *operation;

	if(recursive) {
		iterator = net_getAdjacencyComponentIterator(net);
		while((adjacencyComponent = net_getNextAdjacencyComponent(iterator)) != NULL) {
			net_destruct(adjacencyComponent_getNestedNet(adjacencyComponent), recursive);
		}
		net_destructAdjacencyComponentIterator(iterator);
	}

	netDisk_unloadNet(net->netDisk, net);

	avl_destroy(net->sequences, NULL);

	while((end = net_getFirstEnd(net)) != NULL) {
		end_destruct(end);
	}
	avl_destroy(net->ends, NULL);

	while((atom = net_getFirstAtom(net)) != NULL) {
		atom_destruct(atom);
	}
	avl_destroy(net->atoms, NULL);

	while((adjacencyComponent = net_getFirstAdjacencyComponent(net)) != NULL) {
		adjacencyComponent_destruct(adjacencyComponent);
	}
	avl_destroy(net->adjacencyComponents, NULL);

	while((chain = net_getFirstChain(net)) != NULL) {
		chain_destruct(chain);
	}
	avl_destroy(net->chains, NULL);

	while((operation = net_getFirstOperation(net)) != NULL) {
		operation_destruct(operation);
	}
	avl_destroy(net->operations, NULL);

	free(net->name);
	if(net->parentNetName != NULL) {
		free(net->parentNetName);
	}
	free(net);
}

void net_addSequence(Net *net, Sequence *sequence) {
	avl_insert(net->sequences, sequence);
}

Sequence *net_getFirstSequence(Net *net) {
	return getFirst(net->sequences);
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
	return getIterator(net->sequences);
}

Sequence *net_getNextSequence(Net_SequenceIterator *sequenceIterator) {
	return getNext(sequenceIterator);
}

void net_destructSequenceIterator(Net_SequenceIterator *sequenceIterator) {
	destructIterator(sequenceIterator);
}

End *net_getFirstEnd(Net *net) {
	return getFirst(net->ends);
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
	return getIterator(net->ends);
}

End *net_getNextEnd(Net_EndIterator *endIterator) {
	return getNext(endIterator);
}

void net_destructEndIterator(Net_EndIterator *endIterator) {
	destructIterator(endIterator);
}

Atom *net_getFirstAtom(Net *net) {
	return getFirst(net->atoms);
}

Atom *net_getAtom(Net *net, const char *name) {
	static Atom atom;
	static struct AtomContents atomContents;
	atom.atomContents = &atomContents;
	atomContents.name = (char *)name;
	return avl_find(net->atoms, &atom);
}

int32_t net_getAtomNumber(Net *net) {
	return avl_count(net->atoms);
}

Net_AtomIterator *net_getAtomIterator(Net *net) {
	return getIterator(net->atoms);
}

Atom *net_getNextAtom(Net_AtomIterator *atomIterator) {
	return getNext(atomIterator);
}

void net_destructAtomIterator(Net_AtomIterator *atomIterator) {
	destructIterator(atomIterator);
}

AdjacencyComponent *net_getFirstAdjacencyComponent(Net *net) {
	return getFirst(net->adjacencyComponents);
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
	return getIterator(net->adjacencyComponents);
}

AdjacencyComponent *net_getNextAdjacencyComponent(Net_AdjacencyComponentIterator *adjacencyComponentIterator) {
	return getNext(adjacencyComponentIterator);
}

void net_destructAdjacencyComponentIterator(Net_AdjacencyComponentIterator *adjacencyComponentIterator) {
	destructIterator(adjacencyComponentIterator);
}

AdjacencyComponent *net_getParentAdjacencyComponent(Net *net) {
	Net *net2;
	net2 = netDisk_getNet(net_getNetDisk(net), net->parentNetName);
	return net_getAdjacencyComponent(net2, net_getName(net));
}

Chain *net_getFirstChain(Net *net) {
	return getFirst(net->chains);
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
	return getIterator(net->chains);
}

Chain *net_getNextChain(Net_ChainIterator *chainIterator) {
	return getNext(chainIterator);
}

void net_destructChainIterator(Net_ChainIterator *chainIterator) {
	destructIterator(chainIterator);
}

Operation *net_getFirstOperation(Net *net) {
	return getFirst(net->operations);
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
	return getIterator(net->operations);
}

Operation *net_getNextOperation(Net_OperationIterator *operationIterator) {
	return getNext(operationIterator);
}

void net_destructOperationIterator(Net_OperationIterator *operationIterator) {
	destructIterator(operationIterator);
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

void netDisk_addSequence(NetDisk *netDisk, Sequence *sequence) {
}


void netDisk_unloadSequence(NetDisk *netDisk, Sequence *sequence) {
}

int32_t netDisk_netIsLoaded(NetDisk *netDisk, Net *net) {
}

void netDisk_addNet(NetDisk *netDisk, Net *net) {
}

void netDisk_unloadNet(NetDisk *netDisk, Net *net) {
}

