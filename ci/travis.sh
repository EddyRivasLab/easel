#!/usr/bin/env bash

set -e

TMP_DIR=$(mktemp -d)
HOME_BIN=$HOME/bin

sandbox_run_setup()
{
  test -d $HOME_BIN || mkdir $HOME_BIN

  curl https://raw.githubusercontent.com/horta/port-of-hmmer/master/ci/sandbox_run \
      --output $HOME_BIN/sandbox_run
  chmod +x $HOME_BIN/sandbox_run

  hash -r
}

ppc_run()
{
  sshpass -p "root" ssh -t -oStrictHostKeyChecking=no 127.0.0.1 -p 22125 -l root "$@"
}

ppc_start()
{
  PPC_FILE=$1
  NOHUP_OUT=$TMP_DIR/nohup.out
  touch $NOHUP_OUT
  DIR=$(realpath $TRAVIS_BUILD_DIR/../)
  OPTS="-nographic -vga none -L bios \
    -hda $PPC_FILE -m 512M -net user,hostfwd=tcp::22125-:22 -net nic \
    -virtfs local,path=$DIR,mount_tag=host0,security_model=passthrough,id=host0"

  # Start the emulator
  nohup qemu-system-ppc $OPTS >$NOHUP_OUT 2>&1 &
}

ppc_user_setup()
{
  usr_id=$(id -u)
  grp_id=$(id -g)
  usr_name=$(id -un)

  ppc_run groupadd -g $grp_id $usr_name
  ppc_run useradd -u $usr_id -g $grp_id -m $usr_name
  ppc_run "echo $usr_name:$usr_name | chpasswd"
}

ppc_setup()
{
  echo "Setting PPC up..."

  PPC_ZFILE=debian-wheezy-powerpc.qcow2.bz2
  curl http://rest.s3for.me/hmmer/$PPC_ZFILE --output $TMP_DIR/$PPC_ZFILE
  bunzip2 -v $TMP_DIR/$PPC_ZFILE

  ppc_start $TMP_DIR/debian-wheezy-powerpc.qcow2

  # Wait for PPC machine start up to finish.
  tail -f $TMP_DIR/nohup.out | tee /dev/tty | while read LOGLINE
  do
    [[ "${LOGLINE}" == *"Debian GNU/Linux 7 debian-powerpc"* ]] && break
  done

  ppc_user_setup
  opts="umask=0,trans=virtio,version=9p2000.L"
  ppc_run "mount -t 9p -o $opts host0 /hostdir"

  echo "PPC setup is done."
}

sandbox_run_setup

if [ "$TARGET" == "powerpc-unknown-linux-gnu" ]; then
  ppc_setup
fi
