# Simulador del Dogma Central de la Biología Molecular

Aplicación desarrollada en Python para simular de forma interactiva los principales procesos asociados al dogma central de la biología molecular:

**ADN → ADN → ARN → proteína**

El programa permite visualizar las etapas de **replicación, transcripción y traducción**, utilizando **Biopython** para el tratamiento de las secuencias biológicas y **CustomTkinter** para la interfaz gráfica.

---

## Objetivo

El objetivo del proyecto es representar de forma clara e integrada el flujo de la información genética desde una molécula de ADN hasta la síntesis de una proteína.

El simulador permite observar:

- La estructura de la doble cadena de ADN.
- La replicación del ADN.
- La formación de la cadena líder.
- La formación de la cadena rezagada.
- Los fragmentos de Okazaki.
- La participación de las principales enzimas.
- La transcripción del ADN a ARN mensajero.
- La traducción del ARNm a una secuencia de aminoácidos.
- La detección de los codones de inicio y terminación.

---

## Tecnologías utilizadas

- **Python**
- **Biopython**
- **CustomTkinter**
- **Tkinter**

Biopython se utiliza para trabajar con las secuencias de ADN, ARN y proteínas mediante operaciones como:

```python
complement()
reverse_complement()
transcribe()
translate()
```

También permite leer directamente archivos biológicos en formatos como FASTA y GenBank.

---

## Estructura del proyecto

```text
simulador-dogma-central/
│
├── data/
│   ├── lacZ_Ecoli_NC_000913.3.fasta
│   └── lacZ_Ecoli_NC_000913.3.gb
│
├── src/
│   ├── main.py
│   ├── interfaz.py
│   ├── replicacion.py
│   ├── transcripcion.py
│   └── traduccion.py
│
├── requirements.txt
├── .gitignore
└── README.md
```

### `main.py`

Versión del simulador ejecutada por terminal.

Muestra secuencialmente los resultados de replicación, transcripción y traducción.

### `interfaz.py`

Interfaz gráfica del simulador.

Permite navegar de forma interactiva por las distintas etapas del dogma central.

### `replicacion.py`

Contiene la lógica relacionada con la replicación del ADN.

### `transcripcion.py`

Contiene la lógica relacionada con la transcripción del ADN a ARN mensajero.

### `traduccion.py`

Contiene la lógica relacionada con la traducción del ARNm y la síntesis de la proteína.

### `data/`

Contiene las secuencias biológicas utilizadas en la simulación.

---

## Secuencia utilizada

Como secuencia principal se utiliza el gen **lacZ** de *Escherichia coli* K-12 MG1655.

La secuencia contiene:

- **3075 nucleótidos**
- Codón de inicio: `ATG`
- Codón de terminación: `TAA`

Tras la transcripción y traducción se obtiene una proteína de:

**1024 aminoácidos**

Los archivos utilizados se proporcionan en formatos:

```text
FASTA
GenBank
```

---

## Replicación del ADN

Durante la replicación se representa la separación de las dos cadenas de ADN y la síntesis de nuevas hebras.

El simulador muestra la participación de:

- Helicasa.
- Proteínas SSB.
- Primasa.
- ADN polimerasa.
- ADN polimerasa I.
- ADN ligasa.

También permite diferenciar entre:

### Cadena líder

Se sintetiza de manera continua en dirección:

```text
5' → 3'
```

### Cadena rezagada

Se sintetiza de manera discontinua mediante:

```text
Fragmentos de Okazaki
```

Cada fragmento comienza con un pequeño cebador de ARN.

Posteriormente, los fragmentos se unen para formar una cadena continua.

---

## Replicación semiconservativa

Como resultado de la replicación se generan dos moléculas de ADN.

Cada molécula hija contiene:

```text
1 hebra parental
+
1 hebra de nueva síntesis
```

Por este motivo la replicación se denomina **semiconservativa**.

---

## Transcripción

Durante la transcripción se utiliza una de las cadenas del ADN como molde para sintetizar una molécula de ARN mensajero.

La ARN polimerasa lee la cadena molde y genera el ARNm correspondiente.

La principal diferencia entre el ADN y el ARN es la sustitución de:

```text
Timina (T) → Uracilo (U)
```

Por ejemplo:

```text
ADN:
5' ATGACCATGATT 3'

ARNm:
5' AUGACCAUGAUU 3'
```

La transcripción se realiza mediante Biopython:

```python
arn_m = adn.transcribe()
```

---

## Traducción

Durante la traducción, el ribosoma lee la molécula de ARN mensajero en grupos de tres nucleótidos denominados:

```text
codones
```

Por ejemplo:

```text
AUG | ACC | AUG | AUU | ACG
 ↓     ↓     ↓     ↓     ↓
 M     T     M     I     T
```

Cada codón determina un aminoácido.

El simulador también representa la relación entre:

- ARNm.
- Codón.
- ARNt.
- Anticodón.
- Aminoácido.

---

## Codón de inicio

El simulador busca el primer:

```text
AUG
```

que representa el codón de inicio y codifica para:

```text
Metionina (M)
```

---

## Codones de terminación

La traducción finaliza cuando se encuentra uno de los siguientes codones:

```text
UAA
UAG
UGA
```

Estos codones no codifican ningún aminoácido y actúan como señales de terminación.

Para la secuencia `lacZ` utilizada en el proyecto se obtiene:

```text
UAA → STOP
```

---

## Interfaz gráfica

La aplicación incluye una interfaz gráfica desarrollada con CustomTkinter.

Desde ella se puede navegar por:

```text
Inicio
  ↓
ADN
  ↓
Replicación
  ↓
Transcripción
  ↓
Traducción
  ↓
Resultados
```

La interfaz permite:

- Cargar archivos FASTA.
- Cargar archivos GenBank.
- Introducir manualmente una secuencia de ADN.
- Utilizar una secuencia corta de ejemplo.
- Visualizar la doble cadena de ADN.
- Observar la replicación.
- Consultar los fragmentos de Okazaki.
- Visualizar las enzimas implicadas.
- Observar la transcripción.
- Analizar los primeros codones traducidos.
- Consultar los anticodones de los ARNt.
- Visualizar la proteína obtenida.
- Consultar un resumen final de la simulación.

---

## Instalación

### 1. Clonar el repositorio

```bash
git clone https://github.com/BIO-26-27/Practica-1-Simulacion-del-dogma-central-de-la-biologia-molecular.git
```

Entrar en la carpeta:

```bash
cd simulador-dogma-central
```

### 2. Crear un entorno virtual

```bash
python3 -m venv .venv
```

Activarlo en macOS o Linux:

```bash
source .venv/bin/activate
```

En Windows:

```bash
.venv\Scripts\activate
```

### 3. Instalar las dependencias

```bash
pip install -r requirements.txt
```

El archivo `requirements.txt` contiene:

```text
biopython
customtkinter
```

---

## ▶️ Ejecución

Desde la carpeta que contiene el proyecto:

```bash
python simulador-dogma-central/src/interfaz.py
```

Si ya se está dentro de la carpeta:

```bash
python src/interfaz.py
```

También puede ejecutarse la versión por terminal mediante:

```bash
python src/main.py
```

---

## Archivos compatibles

La interfaz permite cargar los siguientes formatos:

### FASTA

```text
.fasta
.fa
.fna
```

### GenBank

```text
.gb
.gbk
.genbank
```

---

## Ejemplo corto

Además de utilizar una secuencia biológica real, la interfaz incluye una secuencia corta para facilitar la explicación del proceso:

```text
ATGAAACCCGGGTTTTAA
```

Separada en codones:

```text
ATG | AAA | CCC | GGG | TTT | TAA
```

Tras la transcripción:

```text
AUG | AAA | CCC | GGG | UUU | UAA
```

Y durante la traducción:

```text
AUG → M
AAA → K
CCC → P
GGG → G
UUU → F
UAA → STOP
```

Este ejemplo permite observar fácilmente el proceso completo sin utilizar una secuencia de miles de nucleótidos.

---

## Simplificaciones del simulador

El programa tiene un objetivo principalmente didáctico, por lo que algunos procesos biológicos se representan de manera simplificada.

Por ejemplo, para visualizar la formación de la cadena rezagada se utilizan fragmentos de Okazaki de:

```text
500 nucleótidos
```

Este tamaño se utiliza únicamente como criterio de simulación y visualización.

La aplicación pretende mostrar de forma clara los conceptos fundamentales del proceso y no reproducir todos los mecanismos moleculares con detalle experimental.

---

## Flujo completo

El programa permite seguir de forma visual el flujo completo:

```text
ADN
 ↓
REPLICACIÓN
 ↓
ADN
 ↓
TRANSCRIPCIÓN
 ↓
ARNm
 ↓
TRADUCCIÓN
 ↓
PROTEÍNA
```

De esta forma se representa de forma integrada el **dogma central de la biología molecular**.
