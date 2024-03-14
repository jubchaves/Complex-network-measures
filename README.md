
Esse arquivo tem por objetivo detalhar as implementações das funções de análise de redes, comprovando sua eficiência. 
1. Degree (grau do vértice)
def graus(hash):
    graus = []
    for i in range(len(hash)):
        graus.append(len(hash[i]))
        
    return graus

Para essa função, utilizou-se um hash, com chave igual ao nome do vértice, que guardava todos os seus vizinhos. Assim, para descobrir a quantos elementos cada vértice se ligava no grafo, bastou utilizar a função len(), dando o número de vizinhos. O output da função em questão é um vetor de tamanho n (sendo n o número de vértices), em que a iésima entrada corresponde ao grau do iésimo vértice. Comparando com a função nx.degree() da biblioteca Networkx, o resultado é o mesmo. 
2. Shortest path length
def spl1(binaria, a, b):
    if a == b:
        return 0
    
    achou = False
    c = 0
    visitados = []
    falta_visitar = [a]
    
    while len(falta_visitar)>0 and not achou:
        linha = falta_visitar.pop()
        c += 1
        for j in range(len(binaria[linha])):
            if binaria[linha][j] == 1:
                if j == b:
                    achou = True
                    break
                elif j not in visitados:
                    falta_visitar.append(j)
                    visitados.append(j)
                
    if not achou:
        return "infinito"
    else:
        return c

Essa função consiste numa busca em largura numa matriz de correlação sem pesos nas arestas. Toma como input, além da matriz, o vértice inicial e o vértice alvo. Se os vértices são o mesmo, retorna zero e se não há caminho entre eles, infinito. Fora isso, a estratégia usada na implementação é percorrer os vizinhos do vértice atual e, caso não sejam o vértice alvo, enfileirá-los para depois percorrer também seus vizinhos segundo o princípio first in first out. A busca é quebrada pelo indicador de passagem “achou”, caso o vértice alvo seja atingido. O output obtido também é coerente com o que se obtém usando a função pronta nx.shortest_path_length(grafo, i, j) da biblioteca Networkx, com a otimização de reconhecer caminhos desconectados, visto que essa última dá erro, ao invés de retornar infinito. 

3. Triângulos 
def triangulos(binaria, i):
    ti = 0
    n = len(binaria[0])
    
    for j in range(n):
        for h in range(n):
            ti += int(binaria[i][j]) * int(binaria[i][h]) *int(binaria[j][h])
    ti = ti/2                 
    return ti

Essa função é uma aplicação direta da fórmula fornecida no artigo de Rubinov. Recebe como input uma matriz sem peso e não direcionada e um vértice i, ao redor do qual serão calculados os triângulos.  Assim, ela usa de dois loops encadeados para verificar se há conexões entre i e cada vizinho j e h e de j e h entre si, caso no qual se qualifica um triângulo. Ao final, esse somatório é dividido por dois para eliminar a ambiguidade do gráfico não direcionado. O output também é o mesmo da função nx.triangles(grafo, i), da biblioteca Networkx. 

 

4. Characteristic Path Length
def cpl(binaria, d, n):
    n = len(d)
    c =0
    for i in range(n):
        for j in range(n):
            try:
                c+= d[(i, j)]
            except: continue
    C = c/n
    return C


Essa função faz, basicamente uma média, dos menores caminhos entre cada par de vértices. Note que para o caso do caminho ser infinito, ela apenas o ignora na média, por meio do bloco try/except.  Para isso ela chama uma função spl(), que é uma expansão da função spl1(), usada para calcular shortest path length, que, ao invés de retornar um inteiro, correspondente ao caminho entre dois vértices específicos, retorna um dicionário com os caminhos entre todos os pares de vértices. Por questões de otimização, optei por não computar caminhos redundantes (i, j) e (j, i) separadamente, de modo que, no dicionário, i <= j, sempre. Note que, por construção, então, basta, fazer uma média simples dos caminhos, dividindo por n, visto que não há redundâncias. O resultado obtido é o mesmo da função nx.average_shortest_path_length(grafo).
def spl(binaria, n):
    caminhos ={}
    v =[]
    for a in range (n):
        for b in range(n):
            if [a,b] in v or [b,a] in v:
                continue
            v.append([a,b])
            if a == b:
                caminhos[a, b] = 0
            else:
                caminhos[a,b] = spl1(binaria, a, b)
    return caminhos
    


5. Eficiência global

def global_efficiency(n, d):
    n = int(n)
    sumi =0
    for i in range(n):
        sumj = 0
        for j in range(n):
            if i==j: continue
            else:
                try:
                    dij =d[i, j]
                except: 
                    dij = d[j, i]
                if dij != "infinito":
                    sumj += 1/dij
        sumi += sumj/(n-1)
    ef = sumi/n
    return ef
                

Essa função consiste basicamente na média das eficiências de cada vértice, sendo essas dadas pelo somatório do inverso do menor caminho entre ele e cada outro  vértice do grafo, dividido por (n-1). Assim, definimos d como o dicionário de menores caminhos obtido via spl() e usamos o bloco try/except, visto que não há redundâncias, portanto, se existe (i, j), não existe (j, i) como chave. Se o caminho é infinito, seu inverso é zero, portanto, não é necessário acrescer nada à soma. Por outro lado, se i == j, ignoramos essa combinação, visto que o somatório definido no artigo, especifica j!= i. Após obtido o somatório em j, encontra-se a média para cada vértice i, dividindo-se por (n-1) e adiciona-se esse resultado ao somatório i. No fim, basta dividir o somatório de i por n para obter a eficiência global. Nesse caso, o valor retornado não é igual ao obtido com a função equivalente da biblioteca Networkx, porém, testei, e eles se igualam caso, no último passo, a divisão seja feita por (n-1) novamente, ao invés de n. Para garantir a coerência com a expressão fornecida no artigo, mantive a divisão por n. 

6. Coeficiente de clustering
#coeficiente de clustering
def clustering(binaria, k):
    c = 0
    for i in range(len(binaria)):
        t = triangulos(binaria, i)
        if k[i] > 1:
            c += (2*t)/(k[i]*(k[i] - 1))
    c = c/len(binaria[0])     
        #else: Ci = 0 => C+=0
    return c

Para o cálculo do coeficiente de clustering, seguindo a formulação do artigo de Rubinov, encontramos, para cada vértice do grafo, o número de triângulos associados a ele, por meio da função triangulos(), previamente descrita, e o seu grau, dado como input por meio do vetor k, obtido a partir da função graus(), em que a iésima entrada de k, corresponde ao grau do vértice i. Essa taxa de clustering para cada vértice é adicionada num somatório, e, por fim, divide-se tudo pelo número de vértices, obtendo, assim, uma média para o grafo. 
7. Transitividade
def transitividade(binaria, k):
    num = 0
    den = 0
    for i in range(len(binaria)):
        ki = k[i]
        t = triangulos(binaria, i)
        num += 2*t
        den += ki*(ki-1)
    try:
        tr = num/den
    except:
        tr = 0
    return tr

A transitividade é calculada segundo a seguinte expressão:

 
Sendo ti o número de triângulos em torno de um vértice i e ki o grau do mesmo vértice i. Para facilitar os cálculos, defini duas variáveis, num (numerador) e den(denominador). Assim, a função toma como input uma matriz de correlação sem peso e não direcionada e um vetor com os graus de cada vértice. Então, chama a função triângulos para cada vértice i e adiciona seu dobro à variável numerador. Calcula também o grau do vértice, por meio do índice na lista de graus e então calcula o denominador como o somatório de ki*(ki-1). Por fim, calcula a transitividade como num/den, num bloco try/except, visto que, se todos os vétices tiverem grau 1, o denominador será zero e, portanto, a divisão ficará indeterminada. Nesse caso, definimos a transitividade como zero. O output é o mesmo da função equivalente na biblioteca Networkx. 
8. Local efficiency 
def localef(b, binaria, k):
    s1 = 0
    for i in range (len(b)):
        new = {}
        vertices = b[i]
        new[i] = vertices
        for a in vertices:
            lista =[i]
            for v in b[a]:
                if v in vertices:
                    lista.append(v)
            new[a] = sorted(lista)
        #caminho só com os vizinhos de i
        d = spl(make_matrix(binaria, new), len(make_matrix(binaria, new)))
        
        for j in range(len(binaria)):
            for h in range(len(binaria)):
                if j == i or j == h:
                    continue
                try: a = d[j, h]
                except: a = d[h, j]
                if a != "infinito" and k[i] >= 2:
                    a = 1/a
                    s1 += binaria[i][j]*binaria[i][h]*a
                else:
                    continue
        if k[i] > 1:
            e2 = s1/(k[i]*(k[i]-1))
        else:
            e2 = 0
        s1 = 0
        
    e3 = e2/(len(binaria[0])) 
    return e3

A eficiência local foi calculada segundo a fórmula:
 
Para calcular djh¬(Ni), ou seja, o menor caminho entre j e h que só tivesse vizinhos de i, a estratégia utilizada foi a criação de grafos menores a partir do hash dado tal que, para cada i, o dado grafo só contivesse i e seus vizinhos, registrando suas respectivas ligações. Assim, após a criação do hash new, ele foi transformado em uma matriz binaria de correlação por meio da função make_matrix, unicamente para possibilitar o uso da função spl() para o cálculo de tais caminhos. Após obter o valor de d dessa maneira, portanto, aplicou-se a fórmula dada, extraindo novamente k do vetor gerado por graus(). 
def make_matrix(binaria, hash):
    M = []
    for i in range(len(binaria)):
        linha= n*[0]
        if i in hash.keys():
            for v in hash[i]:
                linha[v] += 1
        M.append(linha)
    return M

A função make_matrix utilizada, toma uma matriz binaria base e um hash e a partir dele, cria uma nova matriz de correlação com as mesmas dimensões da matriz original, porém apenas com as conexões especificadas no hash. 


9. Modularity 

def modularity(binaria, k, modulos, l):
    Q = 0
    for i in range(len(binaria)):
        for j in range(len(binaria)):
            if modulos[i] == modulos[j]:
                Q += binaria[i][j] - ((k[i]*k[j])/l)
            else: continue
    Q = Q/l
    
    return Q

Para o cálculo da modularidade utilizou-se a seguinte formulação : 
 
Assim, a função recebe como input a matriz de correlação, o vetor com os graus de cada vértice, o vetor com os módulos de cada vértice e o número de links l, definido como 2m, ou seja, o dobro do número de arestas para um gráfico não direcionado. Como delta recebe zero caso mi != mj, portanto, o programa só realiza a somatória caso os módulos sejam iguais. Assim, subtrai-se da conexão aij (1 ou 0) o valor do produto dos graus de cada vértice dividido pelo número de links. Ao final, divide-se o resultado por l novamente. Para testar essa função, não achei análogo em biblioteca pronta, porém, segundo minhas pesquisas, o valor deveria estar entre -1 e 1, o que se provou verdadeiro para os testes realizados. 

10. Closeness centrality 
def closeness(i, d, n):
    c =0
    for j in range(n):
        try: dij = d[i, j]
        except: dij = d[j, i]
        if dij != "infinito":
            c+=dij
    Li = c/(n-1)
    return 1/Li

Closeness centrality é definida como o inverso do caminho característico para um dado vértice i. Assim, seu cálculo é análogo. A função recebe como input um vértice i, o dicionário d com os menores caminhos entre cada par de vértices e o número de vértices n. Assim, para todo j, se ele for diferente de i, achamos o caminho dij no dicionário e, caso ele não seja infinito, o acrescemos a um somatório c. Ao fim do loop, dividimos o resultado por n-1, conforme a fórmula no artigo de Rubinov. Esse valor Li, corresponde ao caminho característico para i, assim, a função retorna 1/Li, que é exatamente closeness centrality. 
11. Betweenness centrality 
def betwenness_centrality(binaria, hash, i, n, d):
    somatorio = 0
    for h in range(len(binaria[0])):
        for j in range(len(binaria[0])):
            if h == j or h == i or j ==i:
                continue
            else:
                try:
                    di = d[j, h]
                except:
                    di = d[h, j]
                if di == 1:
                    phi = 1
                    phii = 0
                else:
                    phi = num_paths(hash, j, h, d)
                    phii = num_paths(hash,i, j,d)*num_paths(hash, i, h, d)
                if phi == 0:
                    continue
                else:
                    somatorio += phii/phi
    r = 1/((n-1)*(n-2))
    bi = r*somatorio
    return bi

A betweennes centrality foi calculada segundo a formula indicada no artigo de Rubinov. Assim, o n foi dado como input na função, bem como a matriz de correlação, o hash correspondente, o vetor d com os menores caminhos e o vértice i sobre o qual ela é calculada. Para calcular o número de menores caminhos phijh utilizou-se a função num_paths, melhor descrita abaixo:

def num_paths(hash, a, b, d):
    visitados =[]
    falta_visitar=[]
    atual = None
    c = 0
    caminhos = {}
    num_caminhos = 0
    
    if a == b:
        return 0
    
    try: di = d[a, b]
    except: di = d[b, a]
    
    falta_visitar.append([a, 0])
    visitados.append([a, 0])
    
    while len(falta_visitar) > 0:
        atual = falta_visitar.pop()
        c = atual[1]
        atual = atual[0]
        c += 1
        
        if c <= di:
            for e in hash[atual]:
                if e == b:
                    num_caminhos +=1
                    caminhos[num_caminhos] = c               
                elif [e, c] not in visitados:
                    falta_visitar.append([e, c])
                    visitados.append([e, c])
        else:
            continue
                        
    return num_caminhos

Essa função conta o número de menores caminhos entre dois pontos por meio de uma busca em profundidade. Desse modo, ela toma o tamanho do menor caminho do vetor d e é por ele limitada, de forma que, se o número de passos de determinado caminho da busca for maior que o caminho mínimo ela para e aquele caminho não é contabilizado, garantindo o isolamento apenas de menores caminhos. Fora isso, é usada a lógica normal de busca em profundidade com a lista falta_visitar() operando como uma pilha da qual são tirados os próximos vértices. No dicionário caminhos, armazena-se todos os caminhos segundo a ordem que foram encontrados e seu tamanho. Embora ele não seja retornado, é útil para, conseguir ver no debug se todos os caminhos foram contabilizados e se são caminhos mínimos. 
Para o cálculo dos caminhos mínimos que passam por i seguimos o princípio simples da combinatória de que eles corresponderão ao produto de caminhos mínimos de j a i pelo número de caminhos mínimos de i a h. Dessa maneira se obtém phi(i) e termina-se a execução da fórmula, retornando a betweennes centrality do nó. 
12. Within-module degree z-score 

def within_module_z_score(hash, n, modulos, i):
    conexoes = {}
    for j in range(n):
        conexoes[j] = []
        modulo = modulos[j]
        for a in hash[j]:
            if modulos[a] == modulo:
                conexoes[j].append(a)
    kmi = len(conexoes[i])
    k = []
    for j in range(n):
        k.append(len(conexoes[j]))
    k_ = statistics.mean(k)
    std =statistics.stdev(k)
    z = (kmi - k_)/std
    
    return z
Tomando como input o vetor de módulos obtido pelo arquivo txt, essa função determina, para cada um dos vizinhos de cada vértice, quais deles estão localizados no mesmo módulo e os lista em “conexões”. Assim, obtém-se um dicionário informando quantas conexões cada vértice tem dentro de seu próprio módulo. Essa distribuição é utilizada para obter essa métrica para o vértice i, mas também para obter a média  e o desvio padrão dessa distribuição, que são calculados a partir das respectivas funções na biblioteca statistics a partir de uma lista que guarda os tamanhos (portanto número de conexões) dos valores associados a cada chave no dicionário. Tendo essas três medidas, obtemos o within module degree z score pela aplicação da fórmula do artigo. 
13. Participation coefficient 
def participation_coefficient(degree, hash, modulos, i):
    soma = 0

    num_mod = len(set(modulos)) #numero de modulos
    for ind in range(num_mod): #para cada modulo
        mi = modulos[ind]
        m =[]
        for a in range(len(modulos)):
            if modulos[a] == mi:
                m.append(a) #achar elementos daquele modulo
        linha = []
        for elemento in hash[i]:
            if elemento in m:
                linha.append(elemento) #elementos do modulo ligados a i
                
        kmi = len(linha)
        ki = degree[i]
        soma += (kmi/ki)**2

    p = 1 - soma  
    
    return p  

Para calcular o coeficiente de participação, a função, primeiramente, coloca a lista de módulos em um set, eliminando repetições e depois calcula len(), obtendo, assim, o número total de módulos distintos existentes no grafo. Agora, cria-se um loop para achar os elementos daque dado módulo. Tendo esses elementos listados, parte-se para os vizinhos de cada vértice, separando aqueles que estão no módulo analisado. Assim, o número de ligações com elementos daquele módulo mi dá kmi, enquanto ki é o grau global daquele vértice. Então, seguindo o artigo, basta adicionar ao somatório o quadrado da razão entre esses dois valores. Após fazer isso para todos os módulos, calculados os coeficientes subtraindo o valor resultante de 1. 
14. Average neighbour degree 
def average_neighbour(hash, i, k):
    soma = 0
    for v in hash[i]:
            soma+= k[v]
    soma = soma/k[i]
    return soma

Esse atributo consiste no grau médio dos vizinhos de um dado vértice i, assim, a função recebe como input um hash cujas chaves são os vértices, e os elementos uma lista dos respectivos vizinhos, o vértice i e o vetor contento os graus de cada vértice k. O cálculo então é feito somando, para cada vizinho de i, o grau daquele vértice e depois dividindo tudo pelo grau de i. Testando a implementação em comparação com a função nx.average_neighbor_degree(grafo) da biblioteca Networkx, notamos que os valores retornados são iguais. 

15. Assortativity coefficient
def assortativity_coefficient(l, k, binaria):
    sum1 = 0
    sum2 = 0
    sum3 = 0
    for i in range(len(binaria[0])):
        for j in range(len(binaria[0])):
            if binaria[i][j] == 1:
                sum1 += k[i]*k[j]
                sum2 += (1/2)*(k[i]+k[j])
                sum3 += (1/2)*(k[i]**2+k[j]**2)
    num = sum1/l - ((sum2/l)**2)
    den = sum3/l - ((sum2/l)**2)
    A = num/den
    return A

O cálculo de assortativity, também é uma aplicação direta da fórmula fornecida no artigo. Os parâmetros tomados pela função são l, o número de links, k, o vetor com os graus de cada vértice e a matriz binária.
Assim, defini três somatórios diferentes: sum1, que computa a soma dos produtos do grau de i pelo grau de j; sum2, que computa a soma da média dos graus dos dois vértices; sum3, que computa a soma da média dos quadrados dos graus dos dois vértices. Depois, apliquei-os na fórmula do artigo conforme abaixo:

 

O que equivale a: 
[sum1/l - ((sum2/l)**2)]/[sum3/l - ((sum2/l)**2)]

Note que separei, no código original, numerador e denominador por pura conveniência, de modo a aumentar a clareza do código. 
Os resultados obtidos batem com os valores retornados pela função equivalente de Networkx.

Nota: Após a escrita desse relatório algumas alterações mínimas foram feitas nas funções de modo a evitar divisões por zero no caso de vértices desconectados. Tais mudanças não alteram o funcionamento geral das funções e podem ser vistas adequadamente no arquivo python. Para as medidas que se tornavam impossíveis para vértices sem vizinhos, as funções retornam zero por padrão.
