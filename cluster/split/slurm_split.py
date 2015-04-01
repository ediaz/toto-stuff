

#l2spl = 'compute[103,105-109,111]'
l2spl = 'compute[049-052,054-061]'

def getnlist(nlist):
  nlist=nlist.split("[")
  if len(nlist) > 1 :
    head=nlist[0]
    tail=nlist[1] 
  else:
    return nlist
  tail=tail.split(",")	
  nlist=[]
  i = 0
  for t in tail:
    tail[i] = t.replace("]","")
    i+=1

  for t in tail:
    elem =  checksplit(head,t)
    nlist += elem
  return nlist

def checksplit(head,element):
  split = element.split("-")
  if len(split)>1:
    start = int(split[0])
    end   = int(split[1])
    return map(lambda x:head+'%d'%x,range(start,end+1))
  else:
    return [head+element]

def getnlist2(nlist):
        nlist=nlist.split("[")
        if len(nlist) > 1 :
                head=nlist[0]
                tail=nlist[1]
        else:
                return nlist
        tail=tail.split(",")
        
        nlist=[]
        for t in tail:
                t=t.replace("]","")
                sub=t.split("-")
                if len(sub) > 1 :
                        start=int(sub[0])
                        end=int(sub[1])+1
                        for s in range(start,end):
                                k="%3.3d" % (s)
                                k=head+k
                                nlist.append(k)
        return nlist

print getnlist(l2spl)
print getnlist2(l2spl)
