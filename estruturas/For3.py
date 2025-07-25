#imprime informações de um estoque

produto = {'nome':'Caneta chic', 'preco':14.99, 'importada':True, 'estoque':600}

for chave in produto:
    print(chave)

for valor in produto.values():
    print(valor)

for chave, valor in produto.items():
    print(chave, '=', valor)