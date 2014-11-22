CC=clang
CFLAGS=-I./Cosmo/include -I./halomodel/include -I./tools/include -I./Coyote/include -I/usr/local/include
LDFLAGS=-L/usr/local/lib -lgsl -lgslcblas -lfftw3

COSMO_DIR=./Cosmo/src
COSMO_OBJS=$(COSMO_DIR)/cosmo.o $(COSMO_DIR)/lensing.o $(COSMO_DIR)/nofz.o $(COSMO_DIR)/reduced_fit.o $(COSMO_DIR)/decomp_eb.o $(COSMO_DIR)/cmb_bao.o $(COSMO_DIR)/lensing_3rd.o $(COSMO_DIR)/sn1a.o

TOOLS_DIR=./tools/src
TOOLS_OBJ=$(TOOLS_DIR)/config.o $(TOOLS_DIR)/errorlist.o $(TOOLS_DIR)/io.o $(TOOLS_DIR)/maths.o $(TOOLS_DIR)/mvdens.o $(TOOLS_DIR)/par.o

HALOMODEL_DIR=./halomodel/src
HALOMODEL_OBJS=$(HALOMODEL_DIR)/halomodel.o $(HALOMODEL_DIR)/hod.o

COYOTE_DIR=./Coyote/src
COYOTE_OBJS=$(COYOTE_DIR)/constants.o $(COYOTE_DIR)/emu.o $(COYOTE_DIR)/hubble.o

##################################################################

test: test.o $(COSMO_OBJS) $(TOOLS_OBJ) $(HALOMODEL_OBJS) $(COYOTE_OBJS)

%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS)

.PHONY: clean
clean:
	rm -f *.o $(COSMO_DIR)/*.o $(TOOLS_DIR)/*.o $(HALOMODEL_DIR)/*.o $(COYOTE_DIR)/*.o test 