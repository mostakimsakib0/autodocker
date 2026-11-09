gcc aclmgr.c -o aclmgr
strip aclmgr

install -Dm755 aclmgr "$PKG/autodocker/aclmgr"
