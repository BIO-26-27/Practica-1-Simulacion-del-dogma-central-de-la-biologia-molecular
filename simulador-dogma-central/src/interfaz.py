import customtkinter as ctk

from tkinter import filedialog
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


# ============================================================
# CONFIGURACIÓN GENERAL
# ============================================================

ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")


class SimuladorDogmaCentral(ctk.CTk):

    def __init__(self):
        super().__init__()

        # ====================================================
        # VENTANA PRINCIPAL
        # ====================================================

        self.title("Simulador del Dogma Central")
        self.geometry("1350x820")
        self.minsize(1100, 680)

        # ====================================================
        # DATOS
        # ====================================================

        self.registro = None
        self.adn = None
        self.ruta_archivo = None

        # Resultados de replicación
        self.cadena_lider = None
        self.cadena_rezagada = None
        self.cebador_lider = None
        self.fragmentos_okazaki = []

        # Resultados de transcripción
        self.arn_m = None

        # Resultados de traducción
        self.proteina = None
        self.codon_stop = None
        self.posicion_stop = None

        # Estado de la simulación
        self.replicacion_realizada = False
        self.transcripcion_realizada = False
        self.traduccion_realizada = False

        # ====================================================
        # GRID PRINCIPAL
        # ====================================================

        self.grid_columnconfigure(0, weight=0)
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        # ====================================================
        # CREAR INTERFAZ
        # ====================================================

        self.crear_barra_lateral()
        self.crear_contenido_principal()

        self.mostrar_inicio()


    # ========================================================
    # BARRA LATERAL
    # ========================================================

    def crear_barra_lateral(self):

        self.sidebar = ctk.CTkFrame(
            self,
            width=220,
            corner_radius=0
        )

        self.sidebar.grid(
            row=0,
            column=0,
            sticky="nsew"
        )

        self.sidebar.grid_propagate(False)
        self.sidebar.pack_propagate(False)

        # ----------------------------------------------------
        # LOGO / TÍTULO
        # ----------------------------------------------------

        ctk.CTkLabel(
            self.sidebar,
            text="🧬 BIOINFO",
            font=ctk.CTkFont(
                size=24,
                weight="bold"
            )
        ).pack(
            padx=20,
            pady=(30, 4)
        )

        ctk.CTkLabel(
            self.sidebar,
            text="Dogma Central",
            text_color="gray70",
            font=ctk.CTkFont(size=14)
        ).pack(
            padx=20,
            pady=(0, 35)
        )

        # ----------------------------------------------------
        # BOTONES
        # ----------------------------------------------------

        self.boton_inicio = self.crear_boton_menu(
            "⌂  Inicio",
            self.mostrar_inicio
        )

        self.boton_adn = self.crear_boton_menu(
            "🧬  ADN",
            self.mostrar_adn
        )

        self.boton_replicacion = self.crear_boton_menu(
            "↻  Replicación",
            self.mostrar_replicacion
        )

        self.boton_transcripcion = self.crear_boton_menu(
            "→  Transcripción",
            self.mostrar_transcripcion
        )

        self.boton_traduccion = self.crear_boton_menu(
            "⚗  Traducción",
            self.mostrar_traduccion
        )

        self.boton_resultados = self.crear_boton_menu(
            "✓  Resultados",
            self.mostrar_resultados
        )

        self.botones_menu = [
            self.boton_inicio,
            self.boton_adn,
            self.boton_replicacion,
            self.boton_transcripcion,
            self.boton_traduccion,
            self.boton_resultados
        ]

        # ----------------------------------------------------
        # TEXTO INFERIOR
        # ----------------------------------------------------

        ctk.CTkLabel(
            self.sidebar,
            text="Práctica 1\nBioinformática",
            justify="left",
            text_color="gray60",
            font=ctk.CTkFont(size=12)
        ).pack(
            side="bottom",
            padx=20,
            pady=25,
            anchor="w"
        )


    def crear_boton_menu(self, texto, comando):

        boton = ctk.CTkButton(
            self.sidebar,
            text=texto,
            command=comando,
            height=46,
            anchor="w",
            fg_color="transparent",
            hover_color=("gray80", "gray25"),
            font=ctk.CTkFont(
                size=15,
                weight="bold"
            )
        )

        boton.pack(
            padx=15,
            pady=4,
            fill="x"
        )

        return boton


    def marcar_activo(self, boton_activo):

        if not hasattr(self, "botones_menu"):
            return

        for boton in self.botones_menu:
            boton.configure(
                fg_color="transparent"
            )

        boton_activo.configure(
            fg_color=("gray75", "gray25")
        )


    # ========================================================
    # CONTENIDO PRINCIPAL
    # ========================================================

    def crear_contenido_principal(self):

        self.main = ctk.CTkFrame(
            self,
            corner_radius=0,
            fg_color=("gray94", "gray10")
        )

        self.main.grid(
            row=0,
            column=1,
            sticky="nsew"
        )

        self.main.grid_columnconfigure(0, weight=1)
        self.main.grid_rowconfigure(0, weight=1)


    def limpiar_pantalla(self):

        for widget in self.main.winfo_children():
            widget.destroy()


    # ========================================================
    # UTILIDADES VISUALES
    # ========================================================

    def crear_encabezado(
        self,
        padre,
        titulo,
        subtitulo
    ):

        ctk.CTkLabel(
            padre,
            text=titulo,
            font=ctk.CTkFont(
                size=32,
                weight="bold"
            )
        ).pack(anchor="w")

        ctk.CTkLabel(
            padre,
            text=subtitulo,
            text_color="gray65",
            font=ctk.CTkFont(size=15)
        ).pack(
            anchor="w",
            pady=(5, 25)
        )


    def crear_caja_secuencia(
        self,
        padre,
        texto,
        altura=150
    ):

        caja = ctk.CTkTextbox(
            padre,
            height=altura,
            font=("Courier", 13),
            wrap="none"
        )

        caja.pack(
            fill="x",
            padx=25,
            pady=(5, 22)
        )

        caja.insert(
            "1.0",
            texto
        )

        caja.configure(
            state="disabled"
        )

        return caja


    # ========================================================
    # INICIO
    # ========================================================

    def mostrar_inicio(self):

        self.marcar_activo(self.boton_inicio)
        self.limpiar_pantalla()

        contenedor = ctk.CTkScrollableFrame(
            self.main,
            fg_color="transparent"
        )

        contenedor.pack(
            fill="both",
            expand=True,
            padx=40,
            pady=35
        )

        self.crear_encabezado(
            contenedor,
            "Simulador del Dogma Central",
            (
                "Simulación interactiva de replicación, "
                "transcripción y traducción mediante Biopython."
            )
        )

        # ----------------------------------------------------
        # CARGA DEL ARCHIVO
        # ----------------------------------------------------

        carga = ctk.CTkFrame(
            contenedor,
            corner_radius=15
        )

        carga.pack(
            fill="x",
            pady=(0, 30)
        )

        ctk.CTkLabel(
            carga,
            text="Secuencia de entrada",
            font=ctk.CTkFont(
                size=20,
                weight="bold"
            )
        ).pack(
            anchor="w",
            padx=25,
            pady=(22, 5)
        )

        ctk.CTkLabel(
            carga,
            text=self.obtener_nombre_archivo(),
            text_color="gray65",
            font=ctk.CTkFont(size=14)
        ).pack(
            anchor="w",
            padx=25,
            pady=(0, 15)
        )

        ctk.CTkButton(
            carga,
            text="Cargar archivo FASTA / GenBank",
            command=self.cargar_archivo,
            height=42,
            font=ctk.CTkFont(
                size=14,
                weight="bold"
            )
        ).pack(
            anchor="w",
            padx=25,
            pady=(0, 23)
        )

        # ----------------------------------------------------
        # SI HAY ADN
        # ----------------------------------------------------

        if self.adn is not None:

            self.mostrar_tarjetas_resumen(
                contenedor
            )

            # Flujo general
            flujo = ctk.CTkFrame(
                contenedor,
                corner_radius=15
            )

            flujo.pack(
                fill="x",
                pady=(25, 10)
            )

            ctk.CTkLabel(
                flujo,
                text="Flujo de la información genética",
                font=ctk.CTkFont(
                    size=18,
                    weight="bold"
                )
            ).pack(
                anchor="w",
                padx=25,
                pady=(20, 12)
            )

            ctk.CTkLabel(
                flujo,
                text="ADN   →   ADN   →   ARNm   →   PROTEÍNA",
                font=ctk.CTkFont(
                    size=22,
                    weight="bold"
                )
            ).pack(
                anchor="w",
                padx=25,
                pady=(0, 20)
            )

            ctk.CTkButton(
                contenedor,
                text="Comenzar simulación  →",
                command=self.mostrar_adn,
                height=48,
                width=220,
                font=ctk.CTkFont(
                    size=15,
                    weight="bold"
                )
            ).pack(
                anchor="w",
                pady=(20, 10)
            )


    # ========================================================
    # TARJETAS DE RESUMEN
    # ========================================================

    def mostrar_tarjetas_resumen(self, padre):

        tarjetas = ctk.CTkFrame(
            padre,
            fg_color="transparent"
        )

        tarjetas.pack(
            fill="x",
            pady=10
        )

        for columna in range(3):
            tarjetas.grid_columnconfigure(
                columna,
                weight=1
            )

        valor_adn = f"{len(self.adn)} nt"

        if self.transcripcion_realizada:
            valor_arn = f"{len(self.arn_m)} nt"
        else:
            valor_arn = "Pendiente"

        if self.traduccion_realizada:
            valor_proteina = f"{len(self.proteina)} aa"
        else:
            valor_proteina = "Pendiente"

        self.crear_tarjeta(
            tarjetas,
            0,
            "ADN",
            valor_adn,
            "Secuencia cargada"
        )

        self.crear_tarjeta(
            tarjetas,
            1,
            "ARNm",
            valor_arn,
            (
                "Completada"
                if self.transcripcion_realizada
                else "Transcripción"
            )
        )

        self.crear_tarjeta(
            tarjetas,
            2,
            "Proteína",
            valor_proteina,
            (
                "Completada"
                if self.traduccion_realizada
                else "Traducción"
            )
        )


    def crear_tarjeta(
        self,
        padre,
        columna,
        titulo,
        valor,
        descripcion
    ):

        tarjeta = ctk.CTkFrame(
            padre,
            corner_radius=14
        )

        tarjeta.grid(
            row=0,
            column=columna,
            sticky="nsew",
            padx=6
        )

        ctk.CTkLabel(
            tarjeta,
            text=titulo,
            text_color="gray65",
            font=ctk.CTkFont(size=14)
        ).pack(
            anchor="w",
            padx=22,
            pady=(18, 3)
        )

        ctk.CTkLabel(
            tarjeta,
            text=valor,
            font=ctk.CTkFont(
                size=25,
                weight="bold"
            )
        ).pack(
            anchor="w",
            padx=22
        )

        ctk.CTkLabel(
            tarjeta,
            text=descripcion,
            text_color="gray55",
            font=ctk.CTkFont(size=12)
        ).pack(
            anchor="w",
            padx=22,
            pady=(4, 18)
        )


    # ========================================================
    # CARGA DE ARCHIVO
    # ========================================================

    def cargar_archivo(self):

        ruta = filedialog.askopenfilename(
            title="Seleccionar secuencia",
            filetypes=[
                (
                    "Archivos de secuencia",
                    "*.gb *.gbk *.genbank *.fasta *.fa *.fna"
                ),
                (
                    "GenBank",
                    "*.gb *.gbk *.genbank"
                ),
                (
                    "FASTA",
                    "*.fasta *.fa *.fna"
                ),
                (
                    "Todos los archivos",
                    "*.*"
                )
            ]
        )

        if not ruta:
            return

        try:

            extension = Path(ruta).suffix.lower()

            if extension in [
                ".gb",
                ".gbk",
                ".genbank"
            ]:
                formato = "genbank"

            elif extension in [
                ".fasta",
                ".fa",
                ".fna"
            ]:
                formato = "fasta"

            else:
                raise ValueError(
                    "Formato de archivo no compatible."
                )

            registro = SeqIO.read(
                ruta,
                formato
            )

            if len(registro.seq) == 0:
                raise ValueError(
                    "El archivo no contiene una secuencia."
                )

            self.registro = registro
            self.adn = registro.seq.upper()
            self.ruta_archivo = ruta

            self.reiniciar_resultados()

            self.mostrar_inicio()

        except Exception as error:

            self.mostrar_error(
                f"No se pudo cargar el archivo.\n\n{error}"
            )


    def reiniciar_resultados(self):

        self.cadena_lider = None
        self.cadena_rezagada = None
        self.cebador_lider = None
        self.fragmentos_okazaki = []

        self.arn_m = None

        self.proteina = None
        self.codon_stop = None
        self.posicion_stop = None

        self.replicacion_realizada = False
        self.transcripcion_realizada = False
        self.traduccion_realizada = False


    # ========================================================
    # ADN
    # ========================================================

    def mostrar_adn(self):

        if not self.comprobar_adn():
            return

        self.marcar_activo(self.boton_adn)
        self.limpiar_pantalla()

        contenedor = ctk.CTkScrollableFrame(
            self.main,
            fg_color="transparent"
        )

        contenedor.pack(
            fill="both",
            expand=True,
            padx=35,
            pady=30
        )

        self.crear_encabezado(
            contenedor,
            "Molécula de ADN",
            "Información de la secuencia cargada mediante Biopython."
        )

        # ----------------------------------------------------
        # DATOS
        # ----------------------------------------------------

        informacion = ctk.CTkFrame(
            contenedor,
            corner_radius=15
        )

        informacion.pack(
            fill="x",
            pady=(0, 20)
        )

        texto_info = (
            f"ID: {self.registro.id}\n\n"
            f"Descripción: {self.registro.description}\n\n"
            f"Longitud: {len(self.adn)} nucleótidos"
        )

        ctk.CTkLabel(
            informacion,
            text=texto_info,
            justify="left",
            anchor="w",
            wraplength=950,
            font=ctk.CTkFont(size=14)
        ).pack(
            fill="x",
            padx=25,
            pady=22
        )

        # ----------------------------------------------------
        # DOBLE CADENA
        # ----------------------------------------------------

        doble = ctk.CTkFrame(
            contenedor,
            corner_radius=15
        )

        doble.pack(
            fill="x",
            pady=(0, 20)
        )

        ctk.CTkLabel(
            doble,
            text="Doble cadena de ADN",
            font=ctk.CTkFont(
                size=20,
                weight="bold"
            )
        ).pack(
            anchor="w",
            padx=25,
            pady=(20, 5)
        )

        ctk.CTkLabel(
            doble,
            text=(
                "Las dos cadenas son complementarias "
                "y antiparalelas."
            ),
            text_color="gray65",
            font=ctk.CTkFont(size=14)
        ).pack(
            anchor="w",
            padx=25,
            pady=(0, 15)
        )

        longitud_visual = min(
            90,
            len(self.adn)
        )

        codificante = str(
            self.adn[:longitud_visual]
        )

        molde = str(
            self.adn[:longitud_visual].complement()
        )

        texto_adn = (
            "Cadena codificante\n\n"
            f"5'  {codificante}  3'\n"
            f"    {'|' * len(codificante)}\n"
            f"3'  {molde}  5'\n\n"
            "Cadena molde"
        )

        self.crear_caja_secuencia(
            doble,
            texto_adn,
            180
        )

        # ----------------------------------------------------
        # SECUENCIA COMPLETA
        # ----------------------------------------------------

        secuencia = ctk.CTkFrame(
            contenedor,
            corner_radius=15
        )

        secuencia.pack(
            fill="x",
            pady=(0, 20)
        )

        ctk.CTkLabel(
            secuencia,
            text="Secuencia completa",
            font=ctk.CTkFont(
                size=18,
                weight="bold"
            )
        ).pack(
            anchor="w",
            padx=25,
            pady=(18, 5)
        )

        self.crear_caja_secuencia(
            secuencia,
            str(self.adn),
            120
        )

        ctk.CTkButton(
            contenedor,
            text="Continuar a replicación  →",
            command=self.mostrar_replicacion,
            height=45,
            width=220
        ).pack(
            anchor="e",
            pady=(5, 25)
        )


    # ========================================================
    # CÁLCULO DE REPLICACIÓN
    # ========================================================

    def calcular_replicacion(
        self,
        tam_okazaki=500,
        tam_cebador=10
    ):

        if self.replicacion_realizada:
            return

        # Cadena líder
        self.cadena_lider = Seq(
            str(self.adn)
        )

        self.cebador_lider = (
            self.cadena_lider[:tam_cebador]
            .transcribe()
        )

        # Fragmentos molde
        fragmentos_molde = []

        for inicio in range(
            0,
            len(self.adn),
            tam_okazaki
        ):

            fragmentos_molde.append(
                self.adn[
                    inicio:
                    inicio + tam_okazaki
                ]
            )

        self.fragmentos_okazaki = []

        # Se procesan desde el extremo opuesto
        for molde in reversed(
            fragmentos_molde
        ):

            fragmento_nuevo = (
                molde.reverse_complement()
            )

            cebador = (
                fragmento_nuevo[
                    :tam_cebador
                ]
                .transcribe()
            )

            self.fragmentos_okazaki.append(
                {
                    "fragmento":
                        fragmento_nuevo,

                    "cebador":
                        cebador
                }
            )

        cadena_completa = "".join(
            str(datos["fragmento"])
            for datos
            in self.fragmentos_okazaki
        )

        self.cadena_rezagada = Seq(
            cadena_completa
        )

        self.replicacion_realizada = True


    # ========================================================
    # REPLICACIÓN
    # ========================================================

    def mostrar_replicacion(self):

        if not self.comprobar_adn():
            return

        self.calcular_replicacion()

        self.marcar_activo(
            self.boton_replicacion
        )

        self.limpiar_pantalla()

        contenedor = ctk.CTkScrollableFrame(
            self.main,
            fg_color="transparent"
        )

        contenedor.pack(
            fill="both",
            expand=True,
            padx=35,
            pady=30
        )

        self.crear_encabezado(
            contenedor,
            "Replicación del ADN",
            "ADN → ADN · Replicación semiconservativa"
        )

        # ----------------------------------------------------
        # ESQUEMA SUPERIOR
        # ----------------------------------------------------

        proceso = ctk.CTkFrame(
            contenedor,
            corner_radius=15
        )

        proceso.pack(
            fill="x",
            pady=(0, 20)
        )

        ctk.CTkLabel(
            proceso,
            text=(
                "HELICASA   →   PRIMASA   →   "
                "ADN POLIMERASA   →   LIGASA"
            ),
            font=ctk.CTkFont(
                size=17,
                weight="bold"
            )
        ).pack(
            padx=25,
            pady=22
        )

        # ----------------------------------------------------
        # HELICASA
        # ----------------------------------------------------

        tarjeta = ctk.CTkFrame(
            contenedor,
            corner_radius=15
        )

        tarjeta.pack(
            fill="x",
            pady=(0, 18)
        )

        ctk.CTkLabel(
            tarjeta,
            text="1. Apertura de la doble hélice",
            font=ctk.CTkFont(
                size=20,
                weight="bold"
            )
        ).pack(
            anchor="w",
            padx=25,
            pady=(20, 5)
        )

        ctk.CTkLabel(
            tarjeta,
            text=(
                "La helicasa separa las dos cadenas "
                "rompiendo los puentes de hidrógeno."
            ),
            text_color="gray70",
            font=ctk.CTkFont(size=14)
        ).pack(
            anchor="w",
            padx=25,
            pady=(0, 15)
        )

        longitud = min(
            70,
            len(self.adn)
        )

        superior = str(
            self.adn[:longitud]
        )

        inferior = str(
            self.adn[
                :longitud
            ].complement()
        )

        texto = (
            "ADN original:\n\n"
            f"5'  {superior}  3'\n"
            f"    {'|' * len(superior)}\n"
            f"3'  {inferior}  5'\n\n"
            "                         ↓ HELICASA ↓\n\n"
            f"5'  {superior}  3'\n\n"
            f"3'  {inferior}  5'"
        )

        self.crear_caja_secuencia(
            tarjeta,
            texto,
            210
        )

        # ----------------------------------------------------
        # PRIMASA
        # ----------------------------------------------------

        tarjeta = ctk.CTkFrame(
            contenedor,
            corner_radius=15
        )

        tarjeta.pack(
            fill="x",
            pady=(0, 18)
        )

        ctk.CTkLabel(
            tarjeta,
            text="2. Primasa y cebador de ARN",
            font=ctk.CTkFont(
                size=20,
                weight="bold"
            )
        ).pack(
            anchor="w",
            padx=25,
            pady=(20, 5)
        )

        ctk.CTkLabel(
            tarjeta,
            text=(
                "La primasa coloca un cebador de ARN "
                "que permite iniciar la síntesis."
            ),
            text_color="gray70",
            font=ctk.CTkFont(size=14)
        ).pack(
            anchor="w",
            padx=25,
            pady=(0, 12)
        )

        ctk.CTkLabel(
            tarjeta,
            text=(
                f"Cebador ARN:   "
                f"5' {self.cebador_lider} 3'"
            ),
            font=("Courier", 16)
        ).pack(
            anchor="w",
            padx=25,
            pady=(5, 22)
        )

        # ----------------------------------------------------
        # CADENA LÍDER
        # ----------------------------------------------------

        tarjeta = ctk.CTkFrame(
            contenedor,
            corner_radius=15
        )

        tarjeta.pack(
            fill="x",
            pady=(0, 18)
        )

        ctk.CTkLabel(
            tarjeta,
            text="3. Cadena líder",
            font=ctk.CTkFont(
                size=20,
                weight="bold"
            )
        ).pack(
            anchor="w",
            padx=25,
            pady=(20, 5)
        )

        ctk.CTkLabel(
            tarjeta,
            text=(
                "La ADN polimerasa sintetiza la cadena "
                "líder de forma continua en dirección 5' → 3'."
            ),
            text_color="gray70",
            font=ctk.CTkFont(size=14)
        ).pack(
            anchor="w",
            padx=25,
            pady=(0, 15)
        )

        texto_lider = (
            "Cadena molde:\n"
            f"3'  {inferior}  5'\n\n"
            "                     ADN polimerasa  →\n\n"
            "Nueva cadena líder:\n"
            f"5'  {str(self.cadena_lider[:longitud])}  3'"
        )

        self.crear_caja_secuencia(
            tarjeta,
            texto_lider,
            170
        )

        # ----------------------------------------------------
        # CADENA REZAGADA
        # ----------------------------------------------------

        tarjeta = ctk.CTkFrame(
            contenedor,
            corner_radius=15
        )

        tarjeta.pack(
            fill="x",
            pady=(0, 18)
        )

        ctk.CTkLabel(
            tarjeta,
            text="4. Cadena rezagada",
            font=ctk.CTkFont(
                size=20,
                weight="bold"
            )
        ).pack(
            anchor="w",
            padx=25,
            pady=(20, 5)
        )

        ctk.CTkLabel(
            tarjeta,
            text=(
                "La ADN polimerasa solo sintetiza en "
                "dirección 5' → 3'. Por ello, la cadena "
                "rezagada se forma mediante fragmentos "
                "de Okazaki."
            ),
            text_color="gray70",
            justify="left",
            wraplength=950,
            font=ctk.CTkFont(size=14)
        ).pack(
            anchor="w",
            padx=25,
            pady=(0, 15)
        )

        # Contenedor de fragmentos
        fragmentos_frame = ctk.CTkFrame(
            tarjeta,
            fg_color="transparent"
        )

        fragmentos_frame.pack(
            fill="x",
            padx=20,
            pady=(0, 18)
        )

        fragmentos_frame.grid_columnconfigure(
            0,
            weight=1
        )

        fragmentos_frame.grid_columnconfigure(
            1,
            weight=1
        )

        for numero, datos in enumerate(
            self.fragmentos_okazaki,
            start=1
        ):

            fila = (
                (numero - 1) // 2
            )

            columna = (
                (numero - 1) % 2
            )

            fragmento = datos[
                "fragmento"
            ]

            cebador = datos[
                "cebador"
            ]

            caja = ctk.CTkFrame(
                fragmentos_frame,
                corner_radius=10,
                fg_color=(
                    "gray85",
                    "gray20"
                )
            )

            caja.grid(
                row=fila,
                column=columna,
                sticky="nsew",
                padx=6,
                pady=6
            )

            ctk.CTkLabel(
                caja,
                text=(
                    f"Fragmento de Okazaki "
                    f"{numero}"
                ),
                font=ctk.CTkFont(
                    size=14,
                    weight="bold"
                )
            ).pack(
                anchor="w",
                padx=16,
                pady=(12, 3)
            )

            ctk.CTkLabel(
                caja,
                text=f"Cebador: {cebador}",
                text_color="gray65",
                font=("Courier", 11)
            ).pack(
                anchor="w",
                padx=16
            )

            ctk.CTkLabel(
                caja,
                text=(
                    f"5' "
                    f"{str(fragmento[:35])}"
                    f"... 3'"
                ),
                font=("Courier", 11)
            ).pack(
                anchor="w",
                padx=16,
                pady=(3, 2)
            )

            ctk.CTkLabel(
                caja,
                text=(
                    f"{len(fragmento)} "
                    f"nucleótidos"
                ),
                text_color="gray55",
                font=ctk.CTkFont(
                    size=11
                )
            ).pack(
                anchor="w",
                padx=16,
                pady=(0, 12)
            )

        # ----------------------------------------------------
        # LIGASA
        # ----------------------------------------------------

        tarjeta = ctk.CTkFrame(
            contenedor,
            corner_radius=15
        )

        tarjeta.pack(
            fill="x",
            pady=(0, 18)
        )

        ctk.CTkLabel(
            tarjeta,
            text="5. ADN ligasa",
            font=ctk.CTkFont(
                size=20,
                weight="bold"
            )
        ).pack(
            anchor="w",
            padx=25,
            pady=(20, 5)
        )

        ctk.CTkLabel(
            tarjeta,
            text=(
                "La ADN ligasa une los fragmentos "
                "de Okazaki formando una cadena continua."
            ),
            text_color="gray70",
            font=ctk.CTkFont(size=14)
        ).pack(
            anchor="w",
            padx=25,
            pady=(0, 12)
        )

        texto_ligasa = (
            "Cadena rezagada final:\n\n"
            f"5' "
            f"{str(self.cadena_rezagada[:80])}"
            f"... 3'"
        )

        self.crear_caja_secuencia(
            tarjeta,
            texto_ligasa,
            110
        )

        # ----------------------------------------------------
        # RESULTADO
        # ----------------------------------------------------

        resultado = ctk.CTkFrame(
            contenedor,
            corner_radius=15,
            border_width=1
        )

        resultado.pack(
            fill="x",
            pady=(0, 20)
        )

        ctk.CTkLabel(
            resultado,
            text="✓ Replicación completada",
            font=ctk.CTkFont(
                size=20,
                weight="bold"
            )
        ).pack(
            anchor="w",
            padx=25,
            pady=(20, 5)
        )

        ctk.CTkLabel(
            resultado,
            text=(
                "La replicación es semiconservativa: "
                "cada molécula hija conserva una cadena "
                "original y contiene una cadena nueva."
            ),
            wraplength=950,
            justify="left",
            text_color="gray70",
            font=ctk.CTkFont(size=14)
        ).pack(
            anchor="w",
            padx=25,
            pady=(0, 20)
        )

        ctk.CTkButton(
            contenedor,
            text="Continuar a transcripción  →",
            command=self.mostrar_transcripcion,
            height=45,
            width=230
        ).pack(
            anchor="e",
            pady=(5, 25)
        )


    # ========================================================
    # TRANSCRIPCIÓN
    # ========================================================

    def calcular_transcripcion(self):

        if self.transcripcion_realizada:
            return

        self.arn_m = self.adn.transcribe()

        self.transcripcion_realizada = True


    def mostrar_transcripcion(self):

        if not self.comprobar_adn():
            return

        self.calcular_transcripcion()

        self.marcar_activo(
            self.boton_transcripcion
        )

        self.limpiar_pantalla()

        contenedor = ctk.CTkScrollableFrame(
            self.main,
            fg_color="transparent"
        )

        contenedor.pack(
            fill="both",
            expand=True,
            padx=35,
            pady=30
        )

        self.crear_encabezado(
            contenedor,
            "Transcripción",
            "ADN → ARNm · Síntesis de ARN mensajero"
        )

        # ----------------------------------------------------
        # FLUJO
        # ----------------------------------------------------

        flujo = ctk.CTkFrame(
            contenedor,
            corner_radius=15
        )

        flujo.pack(
            fill="x",
            pady=(0, 20)
        )

        ctk.CTkLabel(
            flujo,
            text=(
                "ADN MOLDE   →   ARN POLIMERASA   →   ARNm"
            ),
            font=ctk.CTkFont(
                size=20,
                weight="bold"
            )
        ).pack(
            pady=25
        )

        # ----------------------------------------------------
        # ADN
        # ----------------------------------------------------

        tarjeta = ctk.CTkFrame(
            contenedor,
            corner_radius=15
        )

        tarjeta.pack(
            fill="x",
            pady=(0, 18)
        )

        ctk.CTkLabel(
            tarjeta,
            text="1. Selección de la cadena molde",
            font=ctk.CTkFont(
                size=20,
                weight="bold"
            )
        ).pack(
            anchor="w",
            padx=25,
            pady=(20, 5)
        )

        ctk.CTkLabel(
            tarjeta,
            text=(
                "La ARN polimerasa utiliza una de las "
                "cadenas de ADN como molde."
            ),
            text_color="gray70",
            font=ctk.CTkFont(size=14)
        ).pack(
            anchor="w",
            padx=25,
            pady=(0, 15)
        )

        longitud = min(
            80,
            len(self.adn)
        )

        codificante = str(
            self.adn[:longitud]
        )

        molde = str(
            self.adn[
                :longitud
            ].complement()
        )

        texto_adn = (
            "Cadena codificante:\n"
            f"5'  {codificante}  3'\n\n"
            "Cadena molde:\n"
            f"3'  {molde}  5'"
        )

        self.crear_caja_secuencia(
            tarjeta,
            texto_adn,
            135
        )

        # ----------------------------------------------------
        # ARN POLIMERASA
        # ----------------------------------------------------

        tarjeta = ctk.CTkFrame(
            contenedor,
            corner_radius=15
        )

        tarjeta.pack(
            fill="x",
            pady=(0, 18)
        )

        ctk.CTkLabel(
            tarjeta,
            text="2. ARN polimerasa",
            font=ctk.CTkFont(
                size=20,
                weight="bold"
            )
        ).pack(
            anchor="w",
            padx=25,
            pady=(20, 5)
        )

        ctk.CTkLabel(
            tarjeta,
            text=(
                "La ARN polimerasa lee la cadena molde "
                "y sintetiza una molécula complementaria "
                "de ARN mensajero."
            ),
            text_color="gray70",
            wraplength=950,
            justify="left",
            font=ctk.CTkFont(size=14)
        ).pack(
            anchor="w",
            padx=25,
            pady=(0, 15)
        )

        ctk.CTkLabel(
            tarjeta,
            text="T  →  U",
            font=ctk.CTkFont(
                size=28,
                weight="bold"
            )
        ).pack(
            pady=(5, 20)
        )

        # ----------------------------------------------------
        # ARNm
        # ----------------------------------------------------

        tarjeta = ctk.CTkFrame(
            contenedor,
            corner_radius=15
        )

        tarjeta.pack(
            fill="x",
            pady=(0, 18)
        )

        ctk.CTkLabel(
            tarjeta,
            text="3. ARN mensajero",
            font=ctk.CTkFont(
                size=20,
                weight="bold"
            )
        ).pack(
            anchor="w",
            padx=25,
            pady=(20, 5)
        )

        ctk.CTkLabel(
            tarjeta,
            text=(
                f"Longitud: "
                f"{len(self.arn_m)} nucleótidos"
            ),
            text_color="gray65",
            font=ctk.CTkFont(size=14)
        ).pack(
            anchor="w",
            padx=25,
            pady=(0, 10)
        )

        texto_arn = (
            "Primeros nucleótidos:\n\n"
            f"5' "
            f"{str(self.arn_m[:90])}"
            f"... 3'\n\n"
            "ARNm completo:\n\n"
            f"{self.arn_m}"
        )

        self.crear_caja_secuencia(
            tarjeta,
            texto_arn,
            180
        )

        # ----------------------------------------------------
        # RESULTADO
        # ----------------------------------------------------

        resultado = ctk.CTkFrame(
            contenedor,
            corner_radius=15,
            border_width=1
        )

        resultado.pack(
            fill="x",
            pady=(0, 20)
        )

        ctk.CTkLabel(
            resultado,
            text="✓ Transcripción completada",
            font=ctk.CTkFont(
                size=20,
                weight="bold"
            )
        ).pack(
            anchor="w",
            padx=25,
            pady=(20, 5)
        )

        ctk.CTkLabel(
            resultado,
            text=(
                f"Se ha obtenido un ARNm de "
                f"{len(self.arn_m)} nucleótidos."
            ),
            text_color="gray70",
            font=ctk.CTkFont(size=14)
        ).pack(
            anchor="w",
            padx=25,
            pady=(0, 20)
        )

        ctk.CTkButton(
            contenedor,
            text="Continuar a traducción  →",
            command=self.mostrar_traduccion,
            height=45,
            width=220
        ).pack(
            anchor="e",
            pady=(5, 25)
        )


    # ========================================================
    # TRADUCCIÓN
    # ========================================================

    def calcular_traduccion(self):

        if self.traduccion_realizada:
            return

        self.calcular_transcripcion()

        self.proteina = self.arn_m.translate(
            table=11,
            to_stop=True
        )

        stops = {
            "UAA",
            "UAG",
            "UGA"
        }

        self.codon_stop = None
        self.posicion_stop = None

        for posicion in range(
            0,
            len(self.arn_m) - 2,
            3
        ):

            codon = str(
                self.arn_m[
                    posicion:
                    posicion + 3
                ]
            )

            if codon in stops:

                self.codon_stop = codon

                self.posicion_stop = (
                    posicion // 3
                ) + 1

                break

        self.traduccion_realizada = True


    def mostrar_traduccion(self):

        if not self.comprobar_adn():
            return

        self.calcular_traduccion()

        self.marcar_activo(
            self.boton_traduccion
        )

        self.limpiar_pantalla()

        contenedor = ctk.CTkScrollableFrame(
            self.main,
            fg_color="transparent"
        )

        contenedor.pack(
            fill="both",
            expand=True,
            padx=35,
            pady=30
        )

        self.crear_encabezado(
            contenedor,
            "Traducción y síntesis de proteínas",
            "ARNm → proteína · Lectura de codones"
        )

        # ----------------------------------------------------
        # FLUJO
        # ----------------------------------------------------

        flujo = ctk.CTkFrame(
            contenedor,
            corner_radius=15
        )

        flujo.pack(
            fill="x",
            pady=(0, 20)
        )

        ctk.CTkLabel(
            flujo,
            text=(
                "ARNm   →   RIBOSOMA   →   CODONES   "
                "→   AMINOÁCIDOS   →   PROTEÍNA"
            ),
            font=ctk.CTkFont(
                size=18,
                weight="bold"
            )
        ).pack(
            pady=25
        )

        # ----------------------------------------------------
        # CODONES
        # ----------------------------------------------------

        tarjeta = ctk.CTkFrame(
            contenedor,
            corner_radius=15
        )

        tarjeta.pack(
            fill="x",
            pady=(0, 18)
        )

        ctk.CTkLabel(
            tarjeta,
            text="1. Lectura del ARNm",
            font=ctk.CTkFont(
                size=20,
                weight="bold"
            )
        ).pack(
            anchor="w",
            padx=25,
            pady=(20, 5)
        )

        ctk.CTkLabel(
            tarjeta,
            text=(
                "El ribosoma lee el ARNm de tres "
                "en tres nucleótidos. Cada triplete "
                "se denomina codón."
            ),
            text_color="gray70",
            wraplength=950,
            justify="left",
            font=ctk.CTkFont(size=14)
        ).pack(
            anchor="w",
            padx=25,
            pady=(0, 15)
        )

        codones_frame = ctk.CTkFrame(
            tarjeta,
            fg_color="transparent"
        )

        codones_frame.pack(
            fill="x",
            padx=20,
            pady=(0, 20)
        )

        for columna in range(6):
            codones_frame.grid_columnconfigure(
                columna,
                weight=1
            )

        numero_codones = min(
            18,
            len(self.arn_m) // 3
        )

        for numero in range(
            numero_codones
        ):

            inicio = numero * 3

            codon = self.arn_m[
                inicio:
                inicio + 3
            ]

            aminoacido = codon.translate(
                table=11
            )

            fila = numero // 6
            columna = numero % 6

            caja = ctk.CTkFrame(
                codones_frame,
                corner_radius=10,
                fg_color=(
                    "gray85",
                    "gray20"
                )
            )

            caja.grid(
                row=fila,
                column=columna,
                sticky="nsew",
                padx=5,
                pady=5
            )

            ctk.CTkLabel(
                caja,
                text=str(codon),
                font=ctk.CTkFont(
                    size=16,
                    weight="bold"
                )
            ).pack(
                pady=(12, 2)
            )

            ctk.CTkLabel(
                caja,
                text="↓",
                text_color="gray60"
            ).pack()

            ctk.CTkLabel(
                caja,
                text=str(aminoacido),
                font=ctk.CTkFont(
                    size=18,
                    weight="bold"
                )
            ).pack(
                pady=(2, 12)
            )

        # ----------------------------------------------------
        # PROTEÍNA
        # ----------------------------------------------------

        tarjeta = ctk.CTkFrame(
            contenedor,
            corner_radius=15
        )

        tarjeta.pack(
            fill="x",
            pady=(0, 18)
        )

        ctk.CTkLabel(
            tarjeta,
            text="2. Proteína sintetizada",
            font=ctk.CTkFont(
                size=20,
                weight="bold"
            )
        ).pack(
            anchor="w",
            padx=25,
            pady=(20, 5)
        )

        ctk.CTkLabel(
            tarjeta,
            text=(
                f"Longitud: "
                f"{len(self.proteina)} aminoácidos"
            ),
            text_color="gray65",
            font=ctk.CTkFont(size=14)
        ).pack(
            anchor="w",
            padx=25,
            pady=(0, 10)
        )

        texto_proteina = (
            "Primeros 60 aminoácidos:\n\n"
            f"{self.proteina[:60]}\n\n"
            "Últimos 60 aminoácidos:\n\n"
            f"{self.proteina[-60:]}\n\n"
            "Secuencia completa:\n\n"
            f"{self.proteina}"
        )

        self.crear_caja_secuencia(
            tarjeta,
            texto_proteina,
            220
        )

        # ----------------------------------------------------
        # STOP
        # ----------------------------------------------------

        tarjeta = ctk.CTkFrame(
            contenedor,
            corner_radius=15
        )

        tarjeta.pack(
            fill="x",
            pady=(0, 18)
        )

        ctk.CTkLabel(
            tarjeta,
            text="3. Señal de terminación",
            font=ctk.CTkFont(
                size=20,
                weight="bold"
            )
        ).pack(
            anchor="w",
            padx=25,
            pady=(20, 5)
        )

        if self.codon_stop is not None:

            texto_stop = (
                f"Codón STOP encontrado: "
                f"{self.codon_stop}\n\n"
                f"Posición de codón: "
                f"{self.posicion_stop}\n\n"
                "El ribosoma finaliza la traducción."
            )

        else:

            texto_stop = (
                "No se ha encontrado un codón "
                "STOP en el marco de lectura."
            )

        ctk.CTkLabel(
            tarjeta,
            text=texto_stop,
            justify="left",
            font=ctk.CTkFont(size=15)
        ).pack(
            anchor="w",
            padx=25,
            pady=(5, 22)
        )

        # ----------------------------------------------------
        # RESULTADO
        # ----------------------------------------------------

        resultado = ctk.CTkFrame(
            contenedor,
            corner_radius=15,
            border_width=1
        )

        resultado.pack(
            fill="x",
            pady=(0, 20)
        )

        ctk.CTkLabel(
            resultado,
            text="✓ Traducción completada",
            font=ctk.CTkFont(
                size=20,
                weight="bold"
            )
        ).pack(
            anchor="w",
            padx=25,
            pady=(20, 5)
        )

        ctk.CTkLabel(
            resultado,
            text=(
                f"Se ha sintetizado una proteína de "
                f"{len(self.proteina)} aminoácidos."
            ),
            text_color="gray70",
            font=ctk.CTkFont(size=14)
        ).pack(
            anchor="w",
            padx=25,
            pady=(0, 20)
        )

        ctk.CTkButton(
            contenedor,
            text="Ver resultado final  →",
            command=self.mostrar_resultados,
            height=45,
            width=220
        ).pack(
            anchor="e",
            pady=(5, 25)
        )


    # ========================================================
    # RESULTADOS FINALES
    # ========================================================

    def mostrar_resultados(self):

        if not self.comprobar_adn():
            return

        self.calcular_replicacion()
        self.calcular_transcripcion()
        self.calcular_traduccion()

        self.marcar_activo(
            self.boton_resultados
        )

        self.limpiar_pantalla()

        contenedor = ctk.CTkScrollableFrame(
            self.main,
            fg_color="transparent"
        )

        contenedor.pack(
            fill="both",
            expand=True,
            padx=35,
            pady=30
        )

        self.crear_encabezado(
            contenedor,
            "Resultado de la simulación",
            (
                "Flujo completo de la información genética: "
                "ADN → ADN → ARN → proteína."
            )
        )

        # ----------------------------------------------------
        # TARJETAS
        # ----------------------------------------------------

        tarjetas = ctk.CTkFrame(
            contenedor,
            fg_color="transparent"
        )

        tarjetas.pack(
            fill="x",
            pady=(0, 25)
        )

        for columna in range(3):
            tarjetas.grid_columnconfigure(
                columna,
                weight=1
            )

        self.crear_tarjeta(
            tarjetas,
            0,
            "ADN",
            f"{len(self.adn)} nt",
            "Molécula inicial"
        )

        self.crear_tarjeta(
            tarjetas,
            1,
            "ARNm",
            f"{len(self.arn_m)} nt",
            "Transcripción completada"
        )

        self.crear_tarjeta(
            tarjetas,
            2,
            "Proteína",
            f"{len(self.proteina)} aa",
            "Traducción completada"
        )

        # ----------------------------------------------------
        # FLUJO COMPLETO
        # ----------------------------------------------------

        flujo = ctk.CTkFrame(
            contenedor,
            corner_radius=15
        )

        flujo.pack(
            fill="x",
            pady=(0, 20)
        )

        ctk.CTkLabel(
            flujo,
            text="Flujo completo",
            font=ctk.CTkFont(
                size=20,
                weight="bold"
            )
        ).pack(
            anchor="w",
            padx=25,
            pady=(20, 15)
        )

        ctk.CTkLabel(
            flujo,
            text=(
                "ADN\n"
                "↓\n"
                "REPLICACIÓN\n"
                "↓\n"
                "ADN\n"
                "↓\n"
                "TRANSCRIPCIÓN\n"
                "↓\n"
                "ARNm\n"
                "↓\n"
                "TRADUCCIÓN\n"
                "↓\n"
                "PROTEÍNA"
            ),
            font=ctk.CTkFont(
                size=18,
                weight="bold"
            ),
            justify="center"
        ).pack(
            pady=(0, 25)
        )

        # ----------------------------------------------------
        # DATOS FINALES
        # ----------------------------------------------------

        resumen = ctk.CTkFrame(
            contenedor,
            corner_radius=15
        )

        resumen.pack(
            fill="x",
            pady=(0, 20)
        )

        ctk.CTkLabel(
            resumen,
            text="Resumen",
            font=ctk.CTkFont(
                size=20,
                weight="bold"
            )
        ).pack(
            anchor="w",
            padx=25,
            pady=(20, 12)
        )

        stop_texto = (
            self.codon_stop
            if self.codon_stop
            else "No encontrado"
        )

        texto = (
            f"Secuencia: {self.registro.id}\n\n"
            f"ADN inicial: {len(self.adn)} nucleótidos\n\n"
            f"Fragmentos de Okazaki simulados: "
            f"{len(self.fragmentos_okazaki)}\n\n"
            f"ARNm obtenido: {len(self.arn_m)} nucleótidos\n\n"
            f"Proteína sintetizada: "
            f"{len(self.proteina)} aminoácidos\n\n"
            f"Codón STOP: {stop_texto}"
        )

        ctk.CTkLabel(
            resumen,
            text=texto,
            justify="left",
            font=ctk.CTkFont(size=15)
        ).pack(
            anchor="w",
            padx=25,
            pady=(0, 22)
        )

        # ----------------------------------------------------
        # NOTA
        # ----------------------------------------------------

        nota = ctk.CTkFrame(
            contenedor,
            corner_radius=15
        )

        nota.pack(
            fill="x",
            pady=(0, 20)
        )

        ctk.CTkLabel(
            nota,
            text="Nota sobre la simulación",
            font=ctk.CTkFont(
                size=18,
                weight="bold"
            )
        ).pack(
            anchor="w",
            padx=25,
            pady=(18, 5)
        )

        ctk.CTkLabel(
            nota,
            text=(
                "Para representar visualmente la síntesis "
                "discontinua se utilizan fragmentos de "
                "Okazaki de 500 nucleótidos como "
                "simplificación didáctica."
            ),
            text_color="gray70",
            wraplength=950,
            justify="left",
            font=ctk.CTkFont(size=13)
        ).pack(
            anchor="w",
            padx=25,
            pady=(0, 20)
        )

        # ----------------------------------------------------
        # BOTONES
        # ----------------------------------------------------

        botones = ctk.CTkFrame(
            contenedor,
            fg_color="transparent"
        )

        botones.pack(
            fill="x",
            pady=(5, 30)
        )

        ctk.CTkButton(
            botones,
            text="Volver al inicio",
            command=self.mostrar_inicio,
            height=45,
            width=170
        ).pack(
            side="left",
            padx=(0, 10)
        )

        ctk.CTkButton(
            botones,
            text="Cargar otra secuencia",
            command=self.cargar_archivo,
            height=45,
            width=190
        ).pack(
            side="left"
        )


    # ========================================================
    # UTILIDADES
    # ========================================================

    def comprobar_adn(self):

        if self.adn is None:

            self.mostrar_error(
                "Primero debes cargar una secuencia de ADN."
            )

            return False

        return True


    def mostrar_error(self, mensaje):

        ventana = ctk.CTkToplevel(
            self
        )

        ventana.title("Aviso")
        ventana.geometry("430x220")
        ventana.resizable(
            False,
            False
        )

        ventana.transient(
            self
        )

        ventana.grab_set()

        ctk.CTkLabel(
            ventana,
            text="⚠",
            font=ctk.CTkFont(
                size=35
            )
        ).pack(
            pady=(25, 5)
        )

        ctk.CTkLabel(
            ventana,
            text=mensaje,
            wraplength=360,
            justify="center"
        ).pack(
            padx=25,
            pady=10
        )

        ctk.CTkButton(
            ventana,
            text="Aceptar",
            command=ventana.destroy,
            width=120
        ).pack(
            pady=10
        )


    def obtener_nombre_archivo(self):

        if self.ruta_archivo is None:
            return "Ningún archivo cargado"

        return Path(
            self.ruta_archivo
        ).name


# ============================================================
# EJECUCIÓN
# ============================================================

if __name__ == "__main__":

    app = SimuladorDogmaCentral()
    app.mainloop()