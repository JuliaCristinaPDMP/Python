# Informações dadas pelo enunciado do Prof. Shida
valor = [5, 4, 7, 7]
peso = [5, 6, 8, 4]
capacidadeMochila = 13
quantidade = len(valor)

dp = [[0] * (capacidadeMochila + 1) for _ in range(quantidade + 1)]

for i in range(1, quantidade + 1):
    for w in range(capacidadeMochila + 1):
        if peso[i-1] <= w:
            dp[i][w] = max(dp[i-1][w], valor[i-1] + dp[i-1][w - peso[i-1]])
        else:
            dp[i][w] = dp[i-1][w]

print("Matriz de Programação Dinâmica (DP):")
for linha in dp:
    print(linha)

w = capacidadeMochila
itens_selecionados = []

for i in range(quantidade, 0, -1):
    if dp[i][w] != dp[i-1][w]:
        itens_selecionados.append(i)
        w -= peso[i-1]
        
print(f"\nValor máximo: {dp[quantidade][capacidadeMochila]}")
print(f"Itens selecionados (índices): {sorted(itens_selecionados)}")
