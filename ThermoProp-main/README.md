# ThermoProp

ThermoProp reúne rotinas e exemplos para cálculos termodinâmicos baseados em equações de estado cúbicas. O repositório combina modelos de EoS, leitura de propriedades críticas, cálculo de derivadas (analítico, numérico ou via *automatic differentiation*) e algoritmos de equilíbrio de fases para gerar curvas de ponto de bolha/orvalho e flashes isotérmicos.

## Visão geral
- **Modelos de EoS:** van der Waals, Redlich–Kwong, Soave–Redlich–Kwong e Peng–Robinson, com cálculo de volume molar, energia livre de Gibbs residual, entalpias, capacidades caloríficas e coeficientes de fugacidade e suas derivadas.【F:eos/eos.hpp†L12-L77】
- **Propriedades auxiliares:** constantes SI, tipagens convenientes e integração com o `autodiff` para derivadas automáticas de funções multivariadas.【F:includes.hpp†L13-L38】
- **Banco de dados:** leitura de temperaturas e pressões críticas, além do fator acêntrico, a partir de arquivos YAML simples para compor misturas multicomponentes.【F:database/reading.hpp†L1-L18】
- **Equilíbrio de fases:** implementação dos esquemas de Michelsen/Mollerup para flash isotérmico, rastreio de isocurvas (ponto de bolha, orvalho e fração de vapor especificada) e análise de estabilidade, com opções para jacobiano analítico ou numérico e derivadas de ln φ (analíticas, automáticas ou numéricas).【F:equilibrium/equilibrium.hpp†L1-L85】【F:equilibrium/equilibrium.hpp†L100-L138】
- **Casos de exemplo:** scripts que demonstram a montagem de parâmetros de mistura, a leitura do banco de dados e a geração de curvas de equilíbrio para misturas reais.

### Casos disponíveis em `cases/`
- `analise_derivadas_lnphi/`: calcula derivadas de ln φ por autodiferenciação e por diferenças centrais em relação a T, P e composição para uma mistura de gás natural (C1/C2/iC5), gravando valores e erros em `output.out`. Útil para validar as derivadas analíticas implementadas.【F:cases/analise_derivadas_lnphi/main.cpp†L19-L119】【F:cases/analise_derivadas_lnphi/main.cpp†L136-L191】
- `calculoVolumesMolares/`: exemplo mínimo de carga de banco de dados e cálculo de propriedades para CO2/H2O com Peng–Robinson; prepara `CubicEOSProps` após chamar `compute` e serve como ponto de partida para avaliar volumes molares e fugacidades.【F:cases/calculoVolumesMolares/main.cpp†L1-L24】
- `compare_eos_master/`: compara a saída do código com um *benchmark* externo (EoS Master) para metano puro, reportando fator de compressibilidade e coeficiente de fugacidade em uma condição especificada.【F:cases/compare_eos_master/main.cpp†L6-L38】
- `natural_gas_1/`: traça curvas de ponto de bolha e orvalho para uma mistura típica de gás natural (C1, C2 e iC5), usando diferentes modos de derivadas de ln φ (analítica, automática ou numérica) e registrando o caminho das isocurvas em arquivos `isocurve*.log` e `output_*.out`.【F:cases/natural_gas_1/main.cpp†L1-L53】
- `teste_autodiff/`: valida as derivadas automáticas de ln φ em relação a T, P e composição para um sistema binário CO2/H2O, comparando-as com diferenças finitas e imprimindo o erro relativo no console.【F:cases/teste_autodiff/main.cpp†L1-L67】

## Estrutura do repositório
- `eos/`: funções de EoS cúbica e rotinas auxiliares (cálculo de parâmetros *alpha*, σ, ε, Ω, ψ e estado físico da fase).【F:eos/eos.hpp†L79-L116】
- `equilibrium/`: algoritmos de equilíbrio (flash isotérmico, isocurvas, aproximações inicial/matriciais e opções de iteração de Newton).【F:equilibrium/equilibrium.hpp†L17-L87】【F:equilibrium/equilibrium.hpp†L94-L138】
- `database/`: leitor das propriedades críticas e de composição (`reading.hpp/.cpp`) e arquivos YAML de componentes (`SNG1.yml`, `test.yml`).【F:database/reading.hpp†L1-L18】
- `metodosMatematicos/`: utilidades matemáticas (Newton, manipulação vetorial e integrações numéricas) compartilhadas pelos demais módulos.
- `cases/`: exemplos reproduzíveis. Incluem cálculo de derivadas de ln φ (`analise_derivadas_lnphi`), volumes molares (`calculoVolumesMolares`), comparação de EoS (`compare_eos_master`), traçado de curvas para gás natural (`natural_gas_1`) e testes do `autodiff` (`teste_autodiff`).
- `scripts/`: ferramentas auxiliares, incluindo o *fetch* automático do Eigen.

## Dependências externas
O projeto traz o [`autodiff`](external/autodiff) como dependência interna. Caso o diretório não esteja presente (por exemplo, em um clone raso), basta baixá-lo do repositório oficial:

- Para clonar o `autodiff` diretamente em `external/autodiff`, execute:
  ```bash
  git clone --depth 1 https://github.com/autodiff/autodiff.git external/autodiff
  ```
- Se preferir apontar para outra instalação, ajuste a variável `AutodiffPath` nos `Makefile` dos casos de exemplo para o caminho desejado.

O Eigen pode ser usado a partir da árvore `external` ou do sistema:

- Para baixar e instalar o Eigen localmente em `external/eigen`, execute:
  ```bash
  ./scripts/fetch_eigen.sh
  ```
- Se o Eigen já estiver instalado em outro caminho, defina a variável de ambiente `EIGEN_PATH` ou ajuste o `Makefile` correspondente.
- Caso o diretório `external/eigen` não exista e nenhuma variável seja informada, o `Makefile` utilizará `/usr/include/eigen3`.

## Como estender
- **Novas misturas:** adicione um arquivo YAML em `database/` com `Tc`, `Pc` e `omega` e informe os nomes das espécies ao configurar `read_database`.
- **Modelos de EoS:** selecione qualquer opção do `enum CubicEOSModel` ao montar `param` para testar respostas de diferentes equações cúbicas.【F:eos/eos.hpp†L42-L52】
- **Derivadas:** escolha entre derivadas analíticas, automáticas (`autodiff`) ou numéricas para ln φ por meio de `LnPhiDerivativeType` ao configurar as opções de equilíbrio.【F:equilibrium/equilibrium.hpp†L8-L20】【F:equilibrium/equilibrium.hpp†L55-L86】
- **Algoritmos:** use `computeIsocurve` para traçar bolha/orvalho/fração de vapor ou `flashMichelsen` para flashes isotérmicos com especificações personalizadas.【F:equilibrium/equilibrium.hpp†L104-L126】 Entre em `cases/` para copiar e adaptar exemplos prontos.
