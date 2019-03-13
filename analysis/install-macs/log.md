installing macs (Chen 2009) - done on CentOS server

```bash
cd ~/apps
git clone git@github.com:gchen98/macs.git
cd macs
make all
```

in project folder:

```bash
mkdir bin
cd bin
ln -sv ~/apps/macs/macs .
```
