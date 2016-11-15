containers = cactus ktserver ktremotemgr cpecan

all: ${containers:%=containers.%}

checkclean:
	git diff-index --quiet HEAD


containers.%: checkclean
	cd $* && make build && make push
