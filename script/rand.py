from random import choice

names = [chr(i) for i in range(65,75)]
print "\t".join(["probeset"] + names)
for i in range(10):
	result = [names[i]] + [choice(["0","1"]) for j in range(10)]
	print "\t".join(result)
