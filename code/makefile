objdir=./srcs/
objs=$(objdir)maxima_2d.so $(objdir)flood_shed.so $(objdir)purge_match.so $(objdir)fmm4.so

ifeq ($(MAKECMDGOALS), clean)
	idl_export_dir=. 
endif
ifndef idl_export_dir
	idl_export_dir:=$(shell read -p "Enter direcotry that contains idl_export.h (for example /usr/local/itt/idl81/external/include):" reply; echo $$reply)
endif
cc_command=g++ -fPIC -shared -I$(idl_export_dir) 


all: $(objs)
	@echo "Compiling Done!" 

$(objdir)%.so: $(objdir)%.cpp  
	$(cc_command)$ -o $@ $< 

clean: 
	@rm  $(objdir)*.so ;\
	echo "Clean done!"
