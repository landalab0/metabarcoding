#  :dna: Introducción
El microbioma está formado por comunidades ecológicas de microorganismos que dominan el mundo viviente, en particular, las bacterias las cuales pueden identificarse mediante NGS  (Next Generation Sequencing) con dos vertientes, _Shotgun sequencing:_ donde se secuencia todo el DNA dentro de una muestra para identificar las características funcionales de los genomas y _Metabarcoding:_, donde a través del DNA se revela y cuantifican únicamente la presencia de múltiples especies encontrados en la misma muestra (Taberlet et al., 2018).

También, definimos un _amplicón_ como un fragmento de ADN amplificado mediante la técnica de Reacción en Cadena de la Polimerasa (PCR). Mientras que _metabarcodes_ como una región corta “universal” taxonómicamente informativa flanqueada por dos regiones conservadas que sirven como anclaje para un set de primer universales para ser amplificados por la PCR. Para estos análisis se utilizan genes marcadores como: 

  -> 16S rRNA  bacterias y arqueas 🦠

  -> COI animales 🦈

  ->18S  rRNA  eucariotas 🍄‍🟫

  -> ITS hongos 🍄

  -> rcbl plantas  🌾

La asignación taxonómica de las secuencias generadas es fundamental para la eficacia de estos esfuerzos de biomonitorización molecular. El paso clave en el análisis de secuencias es la manera en que se eliminan los ruidos de los reads y se ensamblan en grupos llamados _ASV_ (Amplicon Sequence Variants) basado en un modelado de errores para secuencias exactas comparado con los _OTU_ (Operational Taxonomic Units) que hacen el clustering con umbral de similitud de ~97%
