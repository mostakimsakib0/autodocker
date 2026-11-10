build()
{
	gcc aclmgr.c -o aclmgr
	strip aclmgr
}

package()
{
	addbin aclmgr autodocker/aclmgr
}
