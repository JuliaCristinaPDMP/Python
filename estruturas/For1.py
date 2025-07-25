#imprime i comecando em 1 ate x<11, vulgo 10
for i in range(1,11):
    print('i = {}'.format(i))

#imprime j 10 vezes com inicio j = 0
for j in range(10):
    print(f'j = {j}')

#imprime a tabuada dos numeros em que x define o numero que está sendo multiplicado
for x in range(1,11):
    for y in range(1,11):
        print(f'{x} * {y} = {x * y}')