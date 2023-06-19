
#
# makefile for NetStart 1.0c
#
# This makefile is used to install the 'netstart' script and manpage
# at proper quarters.
#
# Version: 26 May 2003, K.Rapacki
#

# edit to taste
BINDIR	= /usr/opt/bin
MANDIR	= /usr/opt/man/cat1

# targets, do not change ======================================================
install:	netstart
		cp netstart $(BINDIR)
		chmod 755 $(BINDIR)/netstart

install.man:	man/netstart.1
		cp man/netstart.1 $(MANDIR)
		chmod 644 $(MANDIR)/netstart.1

installc.man:	man/netstart.man
		cat man/netstart.man | compress >$(MANDIR)/netstart.Z
		chmod 644 $(MANDIR)/netstart.Z

# end of makefile =============================================================
