
# coding: utf-8

dic = {}
dic.setdefault('a',[]).append(1) #类似于Perl中的push @{$hash{$k}},$v;
dic.setdefault('a',[]).append(2)
dic.setdefault('b',[]).append(5)
dic.setdefault('b',[]).append(4)

dic

ks =dic.keys()

for k in ks: 
    values= dic[k]
    print(k, values)
    for value in values:
        print (k, value)

