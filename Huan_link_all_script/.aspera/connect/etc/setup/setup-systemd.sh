#!/usr/bin/env sh
#
# Configure Aspera services for use with 'systemd'
#
# Usage: enable | disable
#
# Environment variables:
#   INSTALLDIR - if not present, use /opt/aspera
#
# Aspera, Inc. 2018

_usage()
{
  echo
  echo "  Usage: $0 enable|disable"
  echo
  echo "  Set up Aspera services for use with systemd"
  echo
  exit 1
}

INSTALLDIR=${INSTALLDIR:-/opt/aspera/cargo}
SYSTEMD_LIB=/lib/systemd/system
if [ ! -e $SYSTEMD_LIB ]; then
    SYSTEMD_LIB=/usr/lib/systemd/system
fi

if [ ! -e $SYSTEMD_LIB ]; then
    echo "Can't find systemd lib location"
    exit 1
fi


SYSTEMD_TARGET=/etc/systemd/system/multi-user.target.wants
INITD=/etc/init.d

# Must run as root
if test `id -u` -ne 0; then echo You must be root; exit 1; fi

_initd_clear()
{
    cd $INSTALLDIR/etc/systemd/
    for i in *; do
        rm -f $INITD/${i%.*}
    done
}

_enable()
{
    _initd_clear
    cd $INSTALLDIR/etc/systemd/
    for i in *; do
      cp $i $SYSTEMD_LIB
      ln -sf $SYSTEMD_LIB/$i $SYSTEMD_TARGET
    done
    systemctl daemon-reload
    echo systemd enabled
}

_disable()
{
    cd $INSTALLDIR/etc/systemd/
    for i in *; do
       systemctl stop $i --quiet
       systemctl disable $i  > /dev/null 2>&1
       rm -f $SYSTEMD_LIB/$i
    done
    systemctl daemon-reload
    echo systemd disabled
}

if test $# -ne 1; then
    _usage
fi

if test $1 = enable; then
    _enable
    exit 0
fi

if test $1 = disable; then
    _disable
    exit 0
fi

_usage
exit 1

