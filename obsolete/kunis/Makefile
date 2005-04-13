all:
	make -C lib
	@echo
	@echo "Library completed"
	@echo
	$(foreach dir,$(wildcard example/*),$(MAKE) -C $(dir);)
	@echo
	@echo "Examples completed"
	@echo

install:
	make -C lib
	@echo
	@echo "Library completed"

clean:
	make -C lib clean
	@echo
	@echo "Library clean"
	@echo
	$(foreach dir,$(wildcard example/*),$(MAKE) -C $(dir) clean;)
	@echo
	@echo "Examples clean"
	@echo

rm_makefile:
	rm lib/Makefile
	$(foreach dir,$(wildcard example/*),rm $(dir)/Makefile;)
	@echo
	@echo "Makefiles removed"
	@echo
