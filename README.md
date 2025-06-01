# Regressão beta não linear robusta

Este código foi desenvolvido por ocasião da dissertação de mestrado entitulada **"Regressão beta não linear robusta"** apresentada como requisito parcial para obtenção do título de Mestre em Estatística. 

Com o programa é possível replicar, nas mesmas condições, todos os experimentos, estudos de simulação e aplicações efetuadas por ocasião do referido trabalho, e obter os mesmos resultados.

As funções principais do **[código fonte](https://github.com/eddusousa/nlrobustbetareg/blob/main/Disserta%C3%A7%C3%A3o%20-%20C%C3%B3digos%20e%20Resultados/SourceCode.R)** foram desenvolvidas a partir de adaptações do código base da biblioteca [**robustbetareg**]( https://cran.r-project.org/web/packages/robustbetareg/index.html) do R.

## 🚀 Começando

Não é exigido nenhum processo de instalação, bastando que os arquivos em **[Dissertação - Códigos e Resultados](https://github.com/eddusousa/nlrobustbetareg/tree/main/Disserta%C3%A7%C3%A3o%20-%20C%C3%B3digos%20e%20Resultados)** sejam baixados para o ambiente local para serem abertos e executados por meio do software R, cujas instruções para instalação e uso, além de outras informações, estão disponíveis em https://www.r-project.org/

### 📋 Pré-requisitos

O programa, incluindo as bibliotecas utilizadas, foi desenvolvido e testado utiliando o R versão 4.1.0 e o RStudio versão 2024.12.1 Build 563.

## ⚙️ Executando

### Código Fonte

* As principais funções relacionadas ao ajuste dos modelos e à execução dos experimentos estão no **[código fonte](https://github.com/eddusousa/nlrobustbetareg/blob/main/Disserta%C3%A7%C3%A3o%20-%20C%C3%B3digos%20e%20Resultados/SourceCode.R)**. As principais funções disponíveis são:
  * **robustbetareg**: função principal, responsável por realizar o ajuste do modelo e orquestrar a chamada das demais funções envolvidas na tarefa. Estão disponiveis os ajustes para algumas formas não lineares para os estimadores **MLE**, **SMLE** e **LSMLE**.
  * **FunSample_betareg**: função resposnável por geração de amostra aleatória a partir de um modelo cuja resposta segue uma distribuição beta. 
  * **bootstrap.std.error**: função responsável por calcular o erro padrão via *bootstrap*.
  * **starting.points**: função responsável por obter estimativas iniciais para o processo iterativo de estimação dos parametros.
  * **bootstrap.wald.test**: função resposável por calcular o p-valor via *bootstrap* para o teste tipo-Wald.
  * **select_tuning_constant**: funçãop responsável pelo processo orientado aos dados de seleção da constante de afinação.
  * **apply_spline**: função responsável por aplicar funções pré-definidas para "deslinearizar" a estrutura de regressão.
  * **apply_derivative**: função responsável por derivar a estrutura de regressão segundo a forma não linear escolhida.
  * **simulate_nlrobust**: função responsável por executar o processo de simulação para os cenários pré-definidos.
  * **nlrobust_envelope**: Função responsável por gerar o gráfico de probabilidade normal com envelope simulado. Desenvolvida adaptando a função disponibilizada por **[Terezinha Ribeiro](https://github.com/terezinharibeiro/RobustBetaRegression)**.
  * **run_simdata_experiment**: função responsável por executar o experimento relacionado à aplicação com dados simulados.
  * 

### Simulações

Os códigos utilizados para executar as simulações estão no diretório [Simulações](https://github.com/eddusousa/nlrobustbetareg/tree/main/Disserta%C3%A7%C3%A3o%20-%20C%C3%B3digos%20e%20Resultados/Simula%C3%A7%C3%B5es)

### Aplicações

Os códigos utilizados para executar as aplicações estão no diretório [Aplicações](https://github.com/eddusousa/nlrobustbetareg/tree/main/Disserta%C3%A7%C3%A3o%20-%20C%C3%B3digos%20e%20Resultados/Aplica%C3%A7%C3%B5es)

### Aplicações com dados simulados

[Aplicações com dados simulados](https://github.com/eddusousa/nlrobustbetareg/tree/main/Disserta%C3%A7%C3%A3o%20-%20C%C3%B3digos%20e%20Resultados/Aplica%C3%A7%C3%B5es/Aplica%C3%A7%C3%A3o%20com%20dados%20simulados)

### Aplicações com dados reais

[Aplicações com dados reais](https://github.com/eddusousa/nlrobustbetareg/tree/main/Disserta%C3%A7%C3%A3o%20-%20C%C3%B3digos%20e%20Resultados/Aplica%C3%A7%C3%B5es/Aplica%C3%A7%C3%B5es%20com%20dados%20reais/tuna)

## ✒️ Autor

* **[Eduardo de Sousa Carvalho](https://github.com/eddusousa)**

## 🎁 Colaboração/Apoio

* **[Terezinha Kessia de Assis Ribeiro](https://github.com/terezinharibeiro/)** como orientadora da dissertação de mestrado apoiou todo o desenvolvimento, incluindo a codificação dos experimentos e aplicações, além de fornecer orientações ténicas voltadas para a parte teórica.

