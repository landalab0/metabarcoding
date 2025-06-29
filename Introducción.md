#  :dna: Introducción
El microbioma está formado por comunidades ecológicas de microorganismos que dominan el mundo viviente, en particular, las bacterias las cuales pueden identificarse mediante NGS  (Next Generation Sequencing) con dos vertientes, _Shotgun sequencing:_ donde se secuencia todo el DNA dentro de una muestra para identificar las características funcionales de los genomas, genes funcionales y secuencias no codificantes y _Metabarcoding:_, donde a través del DNA se revela y cuantifican únicamente la presencia de múltiples especies encontrados en la misma muestra.

También, definimos un _amplicón_ como una región específica de ADN o ARN amplificado mediante la técnica de Reacción en Cadena de la Polimerasa (PCR) para estudiar una región de interés. Mientras que un _metabarcode_ es una región corta “universal” taxonómicamente informativa flanqueada por dos secuencias conservadas que sirven como anclaje para un set de primers, y posteriormente se amplifican por la PCR. Para estos análisis se utilizan genes marcadores como: 

  -> 16S rRNA  bacterias y arqueas 🦠

  -> COI animales 🦈

  -> 18S  rRNA eucariotas 🍄‍🟫

  -> ITS hongos 🍄

  -> rcbl plantas  🌾

La asignación taxonómica de las secuencias generadas es fundamental para la eficacia de estos esfuerzos de biomonitorización molecular. Este proceso se facilita mediante el uso de herramientas bioinformáticas para el analisis de secuencias, siendo un paso clave la manera en que se eliminan los ruidos de los reads y se ensamblan en grupos llamados _ASV_ (Amplicon Sequence Variants), basado en un modelado de errores para secuencias exactas que, comparado con los _OTU_ (Operational Taxonomic Units) hacen el clustering con umbral de similitud de ~97%, teniendo los ASV's la ventaja de generar un análisis más detallado en especial en muestras con alta diversidad genética
