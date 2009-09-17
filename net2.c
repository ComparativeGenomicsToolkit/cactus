#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "net2.h"
#include "net2Private.h"

/*
 * Implementation of the net API
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

struct avl_traverser *copyIterator(struct avl_traverser *iterator) {
	struct avl_traverser *copyIterator;
	copyIterator = mallocLocal(sizeof(struct avl_traverser));
	avl_t_copy(copyIterator, iterator);
	return copyIterator;
}

void *getPrevious(struct avl_traverser *iterator) {
	return avl_t_prev(iterator);
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

char *sequence_makeBinaryRepresentation(Sequence *sequence) {

}

char *sequence_makeXMLRepresentation(Sequence *sequence) {

}

Sequence *sequence_loadFromBinaryRepresentation(char *binaryString, NetDisk *netDisk) {

}

Sequence *sequence_loadFromXMLRepresentation(char *xmlString, NetDisk *netDisk) {

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
	//Remove from end.
	end_removeInstance(endInstance_getEnd(endInstance), endInstance);

	destructList(endInstance->children);
	free(endInstance->instance);
	free(endInstance);
}

const char *endInstance_getInstanceName(EndInstance *endInstance) {
	return endInstance->instance;
}

const char *endInstance_getElementName(EndInstance *endInstance) {
	return end_getName(endInstance_getEnd(endInstance));
}

char *endInstance_getCompleteName(EndInstance *endInstance) {
	return netMisc_makeCompleteName(endInstance_getElementName(endInstance), endInstance_getInstanceName(endInstance));
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
	return strcmp(endInstance_getInstanceName((EndInstance *)o1), endInstance_getInstanceName((EndInstance *)o2));
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
			endInstance_constructWithCoordinates(endInstance_getInstanceName(endInstance), end2,
					endInstance_getCoordinate(endInstance), endInstance_getSequence(endInstance));
		}
		else {
			endInstance_construct(endInstance_getInstanceName(endInstance), end2);
		}
	}
	end_destructInstanceIterator(iterator);

	//Copy any parent child links.
	iterator = end_getInstanceIterator(end);
	while((endInstance = end_getNext(iterator)) != NULL) {
		if((endInstance2 = endInstance_getParent(endInstance)) != NULL) {
			endInstance_linkParentAndChild(end_getInstance(end2, endInstance_getInstanceName(endInstance2)),
										   end_getInstance(end2, endInstance_getInstanceName(endInstance)));
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

EndInstance *end_getPrevious(End_InstanceIterator *iterator) {
	return getPrevious(iterator);
}

End_InstanceIterator *end_copyInstanceIterator(End_InstanceIterator *iterator) {
	return copyIterator(iterator);
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
	free(atomInstance->rInstance);
	free(atomInstance);
}

Atom *atomInstance_getAtom(AtomInstance *atomInstance) {
	return atomInstance->atom;
}

const char *atomInstance_getInstanceName(AtomInstance *atomInstance) {
	return endInstance_getInstanceName(atomInstance_getLeft(atomInstance));
}

const char *atomInstance_getElementName(AtomInstance *atomInstance) {
	return atom_getName(atomInstance_getAtom(atomInstance));
}

char *atomInstance_getCompleteName(AtomInstance *atomInstance) {
	return netMisc_makeCompleteName(atomInstance_getInstanceName(atomInstance), atomInstance_getElementName(atomInstance));
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
	return strcmp(atomInstance_getInstanceName((AtomInstance *)o1), atomInstance_getInstanceName((AtomInstance *)o2));
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

AtomInstance *atom_getPrevious(Atom_InstanceIterator *iterator) {
	return getPrevious(iterator);
}

Atom_InstanceIterator *atom_copyInstanceIterator(Atom_InstanceIterator *iterator) {
	return copyIterator(iterator);
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
	//Detach from the parent net.
	net_removeAdjacencyComponent(adjacencyComponent_getNet(adjacencyComponent), adjacencyComponent);
	avl_destroy(adjacencyComponent->ends, NULL);
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

End *adjacencyComponent_getPreviousEnd(AdjacencyComponent_EndIterator *endIterator) {
	return getPrevious(endIterator);
}

AdjacencyComponent_EndIterator *adjacencyComponent_copyEndIterator(AdjacencyComponent_EndIterator *endIterator) {
	return copyIterator(endIterator);
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

Operation *operation_construct(Net *net) {
	Operation *operation;
	operation = malloc(sizeof(Operation));

	net_addOperation(net, operation);
	return operation;
}

void operation_destruct(Operation *operation) {
	net_removeOperation(operation_getNet(operation), operation);
	free(operation);
}

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

Sequence *net_getPreviousSequence(Net_SequenceIterator *sequenceIterator) {
	return getPrevious(sequenceIterator);
}

Net_SequenceIterator *net_copySequenceIterator(Net_SequenceIterator *sequenceIterator) {
	return copyIterator(sequenceIterator);
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

End *net_getPreviousEnd(Net_EndIterator *endIterator) {
	return getPrevious(endIterator);
}

Net_EndIterator *net_copyEndIterator(Net_EndIterator *endIterator) {
	return copyIterator(endIterator);
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

Atom *net_getPreviousAtom(Net_AtomIterator *atomIterator) {
	return getPrevious(atomIterator);
}

Net_AtomIterator *net_copyAtomIterator(Net_AtomIterator *atomIterator) {
	return copyIterator(atomIterator);
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

AdjacencyComponent *net_getPreviousAdjacencyComponent(Net_AdjacencyComponentIterator *adjacencyComponentIterator) {
	return getPrevious(adjacencyComponentIterator);
}

Net_AdjacencyComponentIterator *net_copyAdjacencyComponentIterator(Net_AdjacencyComponentIterator *adjacencyComponentIterator) {
	return copyIterator(adjacencyComponentIterator);
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

Chain *net_getPreviousChain(Net_ChainIterator *chainIterator) {
	return getPrevious(chainIterator);
}

Net_ChainIterator *net_copyChainIterator(Net_ChainIterator *chainIterator) {
	return copyIterator(chainIterator);
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

Operation *net_getPreviousOperation(Net_OperationIterator *operationIterator) {
	return getPrevious(operationIterator);
}

Net_OperationIterator *net_copyOperationIterator(Net_OperationIterator *operationIterator) {
	return copyIterator(operationIterator);
}

void net_destructOperationIterator(Net_OperationIterator *operationIterator) {
	destructIterator(operationIterator);
}

char *net_makeBinaryRepresentation(Net *net) {
}

char *net_makeXMLRepresentation(Net *net) {
}

Net *net_loadFromBinaryRepresentation(char *binaryString, NetDisk *netDisk) {

}

Net *net_loadFromXMLRepresentation(char *xmlString, NetDisk *netDisk) {

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

int32_t netDisk_constructNetsP(const void *o1, const void *o2, void *a) {
	return strcmp(net_getName((Net *)o1), net_getName((Net *)o2));
}

int32_t netDisk_constructSequencesP(const void *o1, const void *o2, void *a) {
	return strcmp(sequence_getName((Sequence *)o1), sequence_getName((Sequence *)o2));
}

NetDisk *netDisk_construct(const char *netDiskFile) {
	NetDisk *netDisk;
	int32_t ecode;

	netDisk = malloc(sizeof(NetDisk));

	//construct lists of in memory objects
	netDisk->sequences = avl_create(net_constructSequencesP, NULL, NULL);
	netDisk->nets = avl_create(net_constructSequencesP, NULL, NULL);

	//make database objects
	netDisk->sequencesDatabase = tcbdbnew();
	netDisk->netsDatabase = tcbdbnew();

	//open the sequences database
	if(!tcbdbopen(netDisk->sequencesDatabase, netDisk->sequencesDatabaseName, BDBOWRITER | BDBOCREAT)){
	   ecode = tcbdbecode(netDisk->sequencesDatabase);
	   fprintf(stderr, "Opening sequences database error: %s\n", tcbdberrmsg(ecode));
	   exit(1);
	}
	//open the nets database
	if(!tcbdbopen(netDisk->netsDatabase, netDisk->netsDatabaseName, BDBOWRITER | BDBOCREAT)){
		ecode = tcbdbecode(netDisk->netsDatabase);
		fprintf(stderr, "Opening nets database error: %s\n", tcbdberrmsg(ecode));
		exit(1);
	}

	return netDisk;
}

void netDisk_destruct(NetDisk *netDisk){
	Sequence *sequence;
	Net *net;
	int32_t ecode;

	while((net = netDisk_getFirstNetInMemory(netDisk)) != NULL) {
		net_destruct(net, FALSE);
	}
	avl_destroy(netDisk->nets, NULL);

	while((sequence = netDisk_getFirstSequenceInMemory(netDisk)) != NULL) {
		sequence_destruct(sequence);
	}
	avl_destroy(netDisk->sequences, NULL);

	//close DBs
	if(!tcbdbclose(netDisk->sequencesDatabase)){
		ecode = tcbdbecode(netDisk->sequencesDatabase);
		fprintf(stderr, "Closing sequences database error: %s\n", tcbdberrmsg(ecode));
		exit(1);
	}
	tcbdbdel(netDisk->sequencesDatabase);

	if(!tcbdbclose(netDisk->netsDatabase)){
		ecode = tcbdbecode(netDisk->netsDatabase);
		fprintf(stderr, "Closing nets database error: %s\n", tcbdberrmsg(ecode));
		exit(1);
	}
	tcbdbdel(netDisk->netsDatabase);

	free(netDisk);
}

int32_t netDisk_write(NetDisk *netDisk){
	NetDisk_NetIterator *netIterator;
	NetDisk_SequenceIterator *sequenceIterator;
	char *cA;
	Net *net;
	Sequence *sequence;
	int32_t ecode;

	netIterator = netDisk_getNetInMemoryIterator(netDisk);
	while((net = netDisk_getNextNet(netIterator)) != NULL) {
		cA = net_makeBinaryRepresentation(net);
		if(!tcbdbput2(netDisk->netsDatabase, net_getName(net), cA)){
			ecode = tcbdbecode(netDisk->netsDatabase);
			fprintf(stderr, "Adding net to database error: %s\n", tcbdberrmsg(ecode));
			return ecode;
		}
		free(cA);
	}
	netDisk_destructNetIterator(netIterator);

	sequenceIterator = netDisk_getSequenceInMemoryIterator(netDisk);
	while((sequence = netDisk_getNextSequence(sequenceIterator)) != NULL) {
		cA = sequence_makeBinaryRepresentation(sequence);
		if(!tcbdbput2(netDisk->sequencesDatabase, net_getName(net), cA)){
			ecode = tcbdbecode(netDisk->sequencesDatabase);
			fprintf(stderr, "Adding sequence to database error: %s\n", tcbdberrmsg(ecode));
			return ecode;
		}
		free(cA);
	}
	netDisk_destructSequenceIterator(sequenceIterator);
	return 0;
}

Sequence *netDisk_getSequence(NetDisk *netDisk, const char *sequenceName) {
	char *cA;
	Sequence *sequence;

	//try in memory list first.
	if((sequence = netDisk_getSequenceInMemory(netDisk, sequenceName)) != NULL) {
		return sequence;
	}
	//else try the database.
	cA = tcbdbget2(netDisk->sequencesDatabase, sequenceName);
	if(cA == NULL) {
		return NULL;
	}
	else {
		sequence = sequence_loadFromBinaryRepresentation(cA, netDisk);
		free(cA);
		return sequence;
	}
}

int32_t netDisk_getSequenceNumberOnDisk(NetDisk *netDisk) {
}

NetDisk_SequenceNameIterator *netDisk_getSequenceNameIterator(NetDisk *netDisk) {
	BDBCUR *iterator;
	iterator = tcbdbcurnew(netDisk->sequencesDatabase);
	tcbdbcurfirst(iterator);
	return iterator;
}

const char *netDisk_getNextSequenceName(NetDisk_SequenceNameIterator *sequenceIterator) {
	return tcbdbcurkey2(sequenceIterator);
}

void netDisk_destructSequenceNameIterator(NetDisk_SequenceNameIterator *sequenceIterator) {
	tcbdbcurdel(sequenceIterator);
}

Sequence *netDisk_getSequenceInMemory(NetDisk *netDisk, const char *sequenceName) {
	static Sequence sequence;
	sequence.name = (char *)sequenceName;
	return avl_find(netDisk->sequences, &sequence);
}

Sequence *netDisk_getFirstSequenceInMemory(NetDisk *netDisk) {
	return getFirst(netDisk->sequences);
}

int32_t netDisk_getSequenceNumberInMemory(NetDisk *netDisk) {
	return avl_count(netDisk->sequences);
}

NetDisk_SequenceIterator *netDisk_getSequenceInMemoryIterator(NetDisk *netDisk) {
	return getIterator(netDisk->sequences);
}

Sequence *netDisk_getNextSequence(NetDisk_SequenceIterator *sequenceIterator) {
	return getNext(sequenceIterator);
}

Sequence *netDisk_getPreviousSequence(NetDisk_SequenceIterator *sequenceIterator) {
	return getPrevious(sequenceIterator);
}

NetDisk_SequenceIterator *netDisk_copySequenceIterator(NetDisk_SequenceIterator *sequenceIterator) {
	return copyIterator(sequenceIterator);
}

void netDisk_destructSequenceIterator(NetDisk_SequenceIterator *sequenceIterator) {
	return destructIterator(sequenceIterator);
}

Net *netDisk_getNet(NetDisk *netDisk, const char *netName) {
}

int32_t netDisk_getNetNumberOnDisk(NetDisk *netDisk) {
}

NetDisk_NetNameIterator *netDisk_getNetNameIterator(NetDisk *netDisk) {
}

const char *netDisk_getNextNetName(NetDisk_NetNameIterator *netIterator) {
}

void netDisk_destructNetNameIterator(NetDisk_NetNameIterator *netIterator) {
}

Net *netDisk_getNetInMemory(NetDisk *netDisk, const char *netName) {
	static Net net;
	net.name = (char *)netName;
	return avl_find(netDisk->nets, &net);
}

Net *netDisk_getFirstNetInMemory(NetDisk *netDisk) {
	return getFirst(netDisk->nets);
}

int32_t netDisk_getNetNumberInMemory(NetDisk *netDisk) {
	return avl_count(netDisk->nets);
}

NetDisk_NetIterator *netDisk_getNetInMemoryIterator(NetDisk *netDisk) {
	return getIterator(netDisk->nets);
}

Net *netDisk_getNextNet(NetDisk_NetIterator *netIterator) {
	return getNext(netIterator);
}

Net *netDisk_getPreviousNet(NetDisk_NetIterator *netIterator) {
	return getPrevious(netIterator);
}

NetDisk_NetIterator *netDisk_copyNetIterator(NetDisk_NetIterator *netIterator) {
	return copyIterator(netIterator);
}

void netDisk_destructNetIterator(NetDisk_NetIterator *netIterator) {
	destructIterator(netIterator);
}

/*
 * Private functions.
 */

int32_t netDisk_deleteSequenceFromDisk(NetDisk *netDisk, const char *sequenceName) {
}

int32_t netDisk_deleteNetFromDisk(NetDisk *netDisk, const char *netName) {
}

void netDisk_addSequence(NetDisk *netDisk, Sequence *sequence) {
#ifdef BEN_DEBUG
	assert(netDisk_getSequence(netDisk) == NULL);
#endif
	avl_insert(netDisk->sequences, sequence);
}

void netDisk_unloadSequence(NetDisk *netDisk, Sequence *sequence) {
#ifdef BEN_DEBUG
	assert(netDisk_getSequenceInMemroy(netDisk) != NULL);
#endif
	avl_delete(netDisk->sequences, sequence);
}

void netDisk_addNet(NetDisk *netDisk, Net *net) {
#ifdef BEN_DEBUG
	assert(netDisk_getNet(netDisk) == NULL);
#endif
	avl_insert(netDisk->nets, net);
}

void netDisk_unloadNet(NetDisk *netDisk, Net *net) {
#ifdef BEN_DEBUG
	assert(netDisk_getNetInMemory(netDisk) == NULL);
#endif
	avl_delete(netDisk->nets, net);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Useful utility functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

char *netMisc_getInstanceName(const char *completeName) {
	char *cA;
	char *cA2;
	char *cA3;
#ifdef BEN_DEBUG
	assert(completeName != NULL);
#endif
	cA = stringCopy(completeName);
	cA2 = cA;
#ifdef BEN_DEBUG
	assert(strlen(cA2) > 0);
#endif
	while(cA2[0] != '.') {
		cA2++;
#ifdef BEN_DEBUG
		assert(strlen(cA2) > 0);
#endif
	}
#ifdef BEN_DEBUG
	assert(cA2[0] == '.');
#endif
	cA3 = stringCopy(cA2+1);
	free(cA);
	return cA3;
}

static char *netMisc_getInstanceNameStatic_instanceName = NULL;
const char *netMisc_getInstanceNameStatic(const char *completeName) {
	if(netMisc_getInstanceNameStatic_instanceName != NULL) {
		free(netMisc_getInstanceNameStatic_instanceName);
	}
	netMisc_getInstanceNameStatic_instanceName = netMisc_getInstanceName(completeName);
	return netMisc_getInstanceNameStatic_instanceName;
}

char *netMisc_getElementName(const char *completeName) {
	char *cA2;
	int32_t i;

	cA2 = (char *)mallocLocal(sizeof(char)*(strlen(completeName)+1));

	i=0;
	while(completeName[i] != '.' && completeName[i] != '\0') {
		cA2[i] = completeName[i];
		i++;
	}

#ifdef BEN_DEBUG
	assert(i <= (int32_t)strlen(completeName));
#endif

	cA2[i] = '\0';
	return cA2;
}

static char *netMisc_getElementNameStatic_elementName = NULL;
const char *netMisc_getElementNameStatic(const char *completeName) {
	if(netMisc_getElementNameStatic_elementName != NULL) {
		free(netMisc_getElementNameStatic_elementName);
	}
	netMisc_getElementNameStatic_elementName = netMisc_getElementName(completeName);
	return netMisc_getElementNameStatic_elementName;
}

char *netMisc_makeCompleteName(const char *elementName, const char *instanceName) {
	char *cA;
	cA = malloc(sizeof(char) * (strlen(elementName) + strlen(instanceName) + 2));
	sprintf(cA, "%s.%s", elementName, instanceName);
	return cA;
}

static char *netMisc_makeCompleteNameStatic_completeName = NULL;
const char *netMisc_makeCompleteNameStatic(const char *elementName, const char *instanceName) {
	if(netMisc_makeCompleteNameStatic_completeName != NULL) {
		free(netMisc_makeCompleteNameStatic_completeName);
	}
	netMisc_makeCompleteNameStatic_completeName = netMisc_makeCompleteName(elementName, instanceName);
	return netMisc_makeCompleteNameStatic_completeName;
}
