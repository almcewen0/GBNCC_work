for ln in open('author_list.txt'):
    l=ln.split('\t')
    nm=l[0]
    if nm[0]=='#':
        continue
    aff=l[1]
    em=l[2]
    orc=l[3].split('\n')[0]
    nmp=''
    for i in range(len(nm.split())):
        if i!=len(nm.split())-1:
            nmp=nmp+nm.split()[i][0]+'.~'
        else:
            nmp=nmp+nm.split()[i]
    if orc=='-':
        nmp='\\author{'+nmp+'}'
    else:
        nmp='\\author['+orc+']'+'{'+nmp+'}'
    print nmp
    for a in aff.split(' AND '):
        print '\\affiliation{'+a+'}'
    print ""


