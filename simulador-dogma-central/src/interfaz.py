import tkinter as tk
from tkinter import filedialog
from pathlib import Path

import customtkinter as ctk
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# ============================================================
# TEMA VISUAL
# ============================================================

ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")

BG = "#111315"
SIDEBAR = "#171A1D"
PANEL = "#1B1F23"
PANEL_2 = "#20252A"
BORDER = "#2D333A"

TEXT = "#F3F4F6"
MUTED = "#A7ADB5"
SUBTLE = "#727A84"

DNA = "#5CA8FF"
RNA = "#58D6B0"
PROTEIN = "#F2C14E"
ACCENT = "#5CA8FF"
SUCCESS = "#58D6B0"
WARNING = "#F2C14E"
ERROR = "#FF7676"

MONO = "Menlo"


AMINOACIDOS = {
    "A": "Alanina",
    "R": "Arginina",
    "N": "Asparagina",
    "D": "Aspartato",
    "C": "Cisteína",
    "E": "Glutamato",
    "Q": "Glutamina",
    "G": "Glicina",
    "H": "Histidina",
    "I": "Isoleucina",
    "L": "Leucina",
    "K": "Lisina",
    "M": "Metionina",
    "F": "Fenilalanina",
    "P": "Prolina",
    "S": "Serina",
    "T": "Treonina",
    "W": "Triptófano",
    "Y": "Tirosina",
    "V": "Valina",
    "*": "STOP",
    "X": "Desconocido",
}


class SimuladorDogmaCentral(ctk.CTk):

    def __init__(self):
        super().__init__()

        self.title("Simulador del Dogma Central")
        self.geometry("1420x860")
        self.minsize(1180, 720)
        self.configure(fg_color=BG)

        # Datos
        self.registro = None
        self.adn = None
        self.ruta_archivo = None

        # Replicación
        self.cadena_lider = None
        self.cadena_rezagada = None
        self.cebador_lider = None
        self.fragmentos_okazaki = []
        self.replicacion_realizada = False

        # Transcripción
        self.arn_m = None
        self.transcripcion_realizada = False

        # Traducción
        self.proteina = None
        self.inicio_traduccion = None
        self.codon_stop = None
        self.posicion_stop_nt = None
        self.posicion_stop_codon = None
        self.codones_traducidos = []
        self.traduccion_realizada = False

        self.grid_columnconfigure(0, weight=0)
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        self.crear_sidebar()
        self.crear_main()
        self.mostrar_inicio()

    # ========================================================
    # ESTRUCTURA GENERAL
    # ========================================================

    def crear_sidebar(self):
        self.sidebar = ctk.CTkFrame(
            self,
            width=215,
            corner_radius=0,
            fg_color=SIDEBAR,
        )
        self.sidebar.grid(row=0, column=0, sticky="nsew")
        self.sidebar.grid_propagate(False)
        self.sidebar.pack_propagate(False)

        marca = ctk.CTkFrame(self.sidebar, fg_color="transparent")
        marca.pack(fill="x", padx=20, pady=(28, 26))

        ctk.CTkLabel(
            marca,
            text="◉",
            text_color=DNA,
            font=ctk.CTkFont(size=25, weight="bold"),
        ).pack(side="left", padx=(0, 9))

        textos = ctk.CTkFrame(marca, fg_color="transparent")
        textos.pack(side="left")

        ctk.CTkLabel(
            textos,
            text="BIOINFO",
            text_color=TEXT,
            font=ctk.CTkFont(size=21, weight="bold"),
        ).pack(anchor="w")

        ctk.CTkLabel(
            textos,
            text="Dogma central",
            text_color=SUBTLE,
            font=ctk.CTkFont(size=12),
        ).pack(anchor="w")

        self.btn_inicio = self.crear_btn_menu("Inicio", self.mostrar_inicio)
        self.btn_adn = self.crear_btn_menu("ADN", self.mostrar_adn)
        self.btn_rep = self.crear_btn_menu("Replicación", self.mostrar_replicacion)
        self.btn_trans = self.crear_btn_menu("Transcripción", self.mostrar_transcripcion)
        self.btn_trad = self.crear_btn_menu("Traducción", self.mostrar_traduccion)
        self.btn_res = self.crear_btn_menu("Resultados", self.mostrar_resultados)

        self.botones_menu = [
            self.btn_inicio,
            self.btn_adn,
            self.btn_rep,
            self.btn_trans,
            self.btn_trad,
            self.btn_res,
        ]

        pie = ctk.CTkFrame(self.sidebar, fg_color="transparent")
        pie.pack(side="bottom", fill="x", padx=20, pady=22)

        ctk.CTkLabel(
            pie,
            text="Práctica 1 · Bioinformática",
            text_color=SUBTLE,
            font=ctk.CTkFont(size=11),
        ).pack(anchor="w")

        ctk.CTkLabel(
            pie,
            text="Python + Biopython",
            text_color=SUBTLE,
            font=ctk.CTkFont(size=11),
        ).pack(anchor="w", pady=(4, 0))

    def crear_btn_menu(self, texto, comando):
        boton = ctk.CTkButton(
            self.sidebar,
            text=texto,
            command=comando,
            height=42,
            anchor="w",
            corner_radius=8,
            fg_color="transparent",
            hover_color=PANEL_2,
            text_color=MUTED,
            font=ctk.CTkFont(size=14, weight="bold"),
        )
        boton.pack(fill="x", padx=12, pady=3)
        return boton

    def marcar_activo(self, activo):
        if not hasattr(self, "botones_menu"):
            return

        for boton in self.botones_menu:
            boton.configure(fg_color="transparent", text_color=MUTED)

        activo.configure(fg_color=PANEL_2, text_color=TEXT)

    def crear_main(self):
        self.main = ctk.CTkFrame(
            self,
            corner_radius=0,
            fg_color=BG,
        )
        self.main.grid(row=0, column=1, sticky="nsew")

    def limpiar(self):
        for widget in self.main.winfo_children():
            widget.destroy()

    def nuevo_contenedor(self, scroll=True):
        if scroll:
            frame = ctk.CTkScrollableFrame(
                self.main,
                fg_color="transparent",
                scrollbar_button_color="#3A4047",
                scrollbar_button_hover_color="#505861",
            )
        else:
            frame = ctk.CTkFrame(self.main, fg_color="transparent")

        frame.pack(fill="both", expand=True, padx=34, pady=26)
        return frame

    # ========================================================
    # COMPONENTES
    # ========================================================

    def encabezado(self, padre, titulo, subtitulo):
        ctk.CTkLabel(
            padre,
            text=titulo,
            text_color=TEXT,
            font=ctk.CTkFont(size=31, weight="bold"),
        ).pack(anchor="w")

        ctk.CTkLabel(
            padre,
            text=subtitulo,
            text_color=MUTED,
            font=ctk.CTkFont(size=14),
        ).pack(anchor="w", pady=(5, 17))

    def panel(self, padre, titulo=None, subtitulo=None, pady=(0, 16)):
        frame = ctk.CTkFrame(
            padre,
            fg_color=PANEL,
            corner_radius=13,
            border_width=1,
            border_color=BORDER,
        )
        frame.pack(fill="x", pady=pady)

        if titulo:
            ctk.CTkLabel(
                frame,
                text=titulo,
                text_color=TEXT,
                font=ctk.CTkFont(size=18, weight="bold"),
            ).pack(anchor="w", padx=22, pady=(18, 2))

        if subtitulo:
            ctk.CTkLabel(
                frame,
                text=subtitulo,
                text_color=MUTED,
                justify="left",
                wraplength=950,
                font=ctk.CTkFont(size=12),
            ).pack(anchor="w", padx=22, pady=(2, 12))

        return frame

    def secuencia_box(self, padre, texto, altura=130):
        caja = ctk.CTkTextbox(
            padre,
            height=altura,
            font=(MONO, 12),
            fg_color="#14171A",
            text_color="#DDE2E7",
            border_width=1,
            border_color=BORDER,
            corner_radius=9,
            wrap="none",
        )
        caja.pack(fill="x", padx=22, pady=(3, 18))
        caja.insert("1.0", texto)
        caja.configure(state="disabled")
        return caja

    def paso_actual(self, padre, etapa):
        nombres = ["ADN", "Replicación", "Transcripción", "Traducción"]

        barra = ctk.CTkFrame(
            padre,
            fg_color=PANEL,
            corner_radius=11,
            border_width=1,
            border_color=BORDER,
        )
        barra.pack(fill="x", pady=(0, 18))

        interior = ctk.CTkFrame(barra, fg_color="transparent")
        interior.pack(pady=12)

        for i, nombre in enumerate(nombres, start=1):
            if i < etapa:
                texto = f"✓  {nombre}"
                color = SUCCESS
            elif i == etapa:
                texto = f"●  {nombre}"
                color = ACCENT
            else:
                texto = f"○  {nombre}"
                color = SUBTLE

            ctk.CTkLabel(
                interior,
                text=texto,
                text_color=color,
                font=ctk.CTkFont(size=13, weight="bold"),
            ).pack(side="left", padx=7)

            if i < len(nombres):
                ctk.CTkLabel(
                    interior,
                    text="›",
                    text_color="#4E555D",
                    font=ctk.CTkFont(size=18),
                ).pack(side="left", padx=2)

    def metrica(self, padre, columna, etiqueta, valor, detalle, color=TEXT):
        tarjeta = ctk.CTkFrame(
            padre,
            fg_color=PANEL,
            corner_radius=12,
            border_width=1,
            border_color=BORDER,
        )
        tarjeta.grid(row=0, column=columna, sticky="nsew", padx=5)

        ctk.CTkLabel(
            tarjeta,
            text=etiqueta,
            text_color=MUTED,
            font=ctk.CTkFont(size=12),
        ).pack(anchor="w", padx=18, pady=(15, 1))

        ctk.CTkLabel(
            tarjeta,
            text=valor,
            text_color=color,
            font=ctk.CTkFont(size=23, weight="bold"),
        ).pack(anchor="w", padx=18)

        ctk.CTkLabel(
            tarjeta,
            text=detalle,
            text_color=SUBTLE,
            font=ctk.CTkFont(size=10),
        ).pack(anchor="w", padx=18, pady=(3, 14))

    def boton_siguiente(self, padre, texto, comando):
        ctk.CTkButton(
            padre,
            text=texto,
            command=comando,
            height=43,
            width=220,
            corner_radius=9,
            fg_color="#246DA8",
            hover_color="#2F80C2",
            font=ctk.CTkFont(size=13, weight="bold"),
        ).pack(anchor="e", pady=(3, 20))

    def formatear(self, secuencia, limite=78):
        s = str(secuencia)
        if len(s) <= limite:
            return s
        return s[:limite] + "…"

    # ========================================================
    # INICIO
    # ========================================================

    def mostrar_inicio(self):
        self.marcar_activo(self.btn_inicio)
        self.limpiar()

        cont = self.nuevo_contenedor()

        self.encabezado(
            cont,
            "Simulador del Dogma Central",
            "Una visualización interactiva del flujo ADN → ARN → proteína.",
        )

        fuente = self.panel(
            cont,
            "Secuencia de entrada",
            "Carga un archivo FASTA o GenBank, utiliza el ejemplo didáctico o introduce una secuencia manual.",
        )

        fila = ctk.CTkFrame(fuente, fg_color="transparent")
        fila.pack(fill="x", padx=22, pady=(0, 12))

        ctk.CTkButton(
            fila,
            text="Cargar FASTA / GenBank",
            command=self.cargar_archivo,
            height=40,
            corner_radius=8,
        ).pack(side="left", padx=(0, 9))

        ctk.CTkButton(
            fila,
            text="Ejemplo corto",
            command=self.cargar_ejemplo_corto,
            height=40,
            corner_radius=8,
            fg_color=PANEL_2,
            hover_color="#30363D",
        ).pack(side="left")

        estado = self.obtener_nombre_archivo()

        ctk.CTkLabel(
            fuente,
            text=estado,
            text_color=SUBTLE,
            font=ctk.CTkFont(size=11),
        ).pack(anchor="w", padx=22, pady=(0, 14))

        self.entrada_manual = ctk.CTkEntry(
            fuente,
            placeholder_text="Secuencia manual, por ejemplo: ATGAAACCCGGGTTTTAA",
            height=40,
            font=(MONO, 12),
            fg_color="#14171A",
            border_color=BORDER,
        )
        self.entrada_manual.pack(fill="x", padx=22, pady=(0, 9))

        ctk.CTkButton(
            fuente,
            text="Usar secuencia manual",
            command=self.cargar_secuencia_manual,
            height=38,
            width=175,
            corner_radius=8,
            fg_color=PANEL_2,
            hover_color="#30363D",
        ).pack(anchor="w", padx=22, pady=(0, 18))

        if self.adn is not None:
            metricas = ctk.CTkFrame(cont, fg_color="transparent")
            metricas.pack(fill="x", pady=(4, 18))

            for i in range(3):
                metricas.grid_columnconfigure(i, weight=1)

            self.metrica(
                metricas,
                0,
                "ADN",
                f"{len(self.adn)} nt",
                "Secuencia cargada",
                DNA,
            )

            arn_valor = f"{len(self.arn_m)} nt" if self.transcripcion_realizada else "—"
            prot_valor = (
                f"{len(self.proteina)} aa"
                if self.traduccion_realizada and self.inicio_traduccion is not None
                else "—"
            )

            self.metrica(metricas, 1, "ARNm", arn_valor, "Tras la transcripción", RNA)
            self.metrica(metricas, 2, "Proteína", prot_valor, "Tras la traducción", PROTEIN)

            resumen = self.panel(cont, "Flujo de trabajo")
            ctk.CTkLabel(
                resumen,
                text="ADN   →   replicación   →   ADN   →   transcripción   →   ARNm   →   traducción   →   proteína",
                text_color=TEXT,
                font=ctk.CTkFont(size=15, weight="bold"),
            ).pack(anchor="w", padx=22, pady=(8, 18))

            ctk.CTkButton(
                cont,
                text="Comenzar simulación",
                command=self.mostrar_adn,
                height=43,
                width=190,
                corner_radius=9,
                font=ctk.CTkFont(size=13, weight="bold"),
            ).pack(anchor="w", pady=(0, 12))

    # ========================================================
    # CARGA
    # ========================================================

    def cargar_archivo(self):
        ruta = filedialog.askopenfilename(
            title="Seleccionar secuencia",
            filetypes=[
                ("Secuencias", "*.gb *.gbk *.genbank *.fasta *.fa *.fna"),
                ("GenBank", "*.gb *.gbk *.genbank"),
                ("FASTA", "*.fasta *.fa *.fna"),
                ("Todos los archivos", "*.*"),
            ],
        )

        if not ruta:
            return

        try:
            ext = Path(ruta).suffix.lower()

            if ext in [".gb", ".gbk", ".genbank"]:
                formato = "genbank"
            elif ext in [".fasta", ".fa", ".fna"]:
                formato = "fasta"
            else:
                raise ValueError("Formato no compatible.")

            registro = SeqIO.read(ruta, formato)

            if len(registro.seq) == 0:
                raise ValueError("La secuencia está vacía.")

            self.registro = registro
            self.adn = registro.seq.upper()
            self.ruta_archivo = ruta

            self.reiniciar_resultados()
            self.mostrar_inicio()

        except Exception as exc:
            self.mostrar_error(f"No se pudo cargar la secuencia.\n\n{exc}")

    def cargar_ejemplo_corto(self):
        seq = Seq("ATGAAACCCGGGTTTTAA")
        self.registro = SeqRecord(
            seq,
            id="EJEMPLO_CORTO",
            description="Secuencia didáctica: inicio AUG y terminación UAA",
        )
        self.adn = seq
        self.ruta_archivo = None
        self.reiniciar_resultados()
        self.mostrar_inicio()

    def cargar_secuencia_manual(self):
        texto = (
            self.entrada_manual.get()
            .upper()
            .replace(" ", "")
            .replace("\n", "")
            .replace("\t", "")
        )

        if not texto:
            self.mostrar_error("Introduce una secuencia de ADN.")
            return

        invalidos = set(texto) - set("ATCG")
        if invalidos:
            self.mostrar_error("La secuencia manual solo puede contener A, T, C y G.")
            return

        self.adn = Seq(texto)
        self.registro = SeqRecord(
            self.adn,
            id="SECUENCIA_MANUAL",
            description="Secuencia introducida manualmente",
        )
        self.ruta_archivo = None
        self.reiniciar_resultados()
        self.mostrar_inicio()

    def reiniciar_resultados(self):
        self.cadena_lider = None
        self.cadena_rezagada = None
        self.cebador_lider = None
        self.fragmentos_okazaki = []
        self.replicacion_realizada = False

        self.arn_m = None
        self.transcripcion_realizada = False

        self.proteina = None
        self.inicio_traduccion = None
        self.codon_stop = None
        self.posicion_stop_nt = None
        self.posicion_stop_codon = None
        self.codones_traducidos = []
        self.traduccion_realizada = False

    # ========================================================
    # ADN
    # ========================================================

    def mostrar_adn(self):
        if not self.comprobar_adn():
            return

        self.marcar_activo(self.btn_adn)
        self.limpiar()
        cont = self.nuevo_contenedor()

        self.encabezado(
            cont,
            "Molécula de ADN",
            "Secuencia de partida y relación entre la hebra codificante y la hebra molde.",
        )
        self.paso_actual(cont, 1)

        datos = self.panel(cont, "Información de la secuencia")

        ctk.CTkLabel(
            datos,
            text=f"ID   {self.registro.id}",
            text_color=TEXT,
            font=ctk.CTkFont(size=13, weight="bold"),
        ).pack(anchor="w", padx=22, pady=(9, 2))

        ctk.CTkLabel(
            datos,
            text=f"{self.registro.description}",
            text_color=MUTED,
            wraplength=1000,
            justify="left",
            font=ctk.CTkFont(size=12),
        ).pack(anchor="w", padx=22, pady=(2, 4))

        ctk.CTkLabel(
            datos,
            text=f"{len(self.adn)} nucleótidos",
            text_color=DNA,
            font=ctk.CTkFont(size=13, weight="bold"),
        ).pack(anchor="w", padx=22, pady=(3, 18))

        doble = self.panel(
            cont,
            "Doble cadena",
            "Las hebras son complementarias y antiparalelas.",
        )

        n = min(72, len(self.adn))
        cod = str(self.adn[:n])
        molde = str(self.adn[:n].complement())

        texto = (
            f"Cadena codificante   5'  {cod}  3'\n"
            f"                         {'│' * len(cod)}\n"
            f"Cadena molde         3'  {molde}  5'\n\n"
            "Complementariedad: A ↔ T    ·    C ↔ G"
        )
        self.secuencia_box(doble, texto, 145)

        nota = self.panel(
            cont,
            "Suposición del simulador",
            "La secuencia cargada se interpreta como la hebra codificante escrita 5' → 3'.",
        )
        ctk.CTkLabel(
            nota,
            text="Esta convención permite obtener la hebra molde mediante complement() y el ARNm mediante transcribe().",
            text_color=MUTED,
            font=ctk.CTkFont(size=12),
        ).pack(anchor="w", padx=22, pady=(2, 18))

        self.boton_siguiente(cont, "Continuar a replicación  →", self.mostrar_replicacion)

    # ========================================================
    # REPLICACIÓN
    # ========================================================

    def calcular_replicacion(self, tam_okazaki=500, tam_cebador=10):
        if self.replicacion_realizada:
            return

        self.cadena_lider = Seq(str(self.adn))

        n_cebador = min(tam_cebador, len(self.adn))
        self.cebador_lider = self.adn[:n_cebador].transcribe()

        fragmentos_molde = [
            self.adn[i:i + tam_okazaki]
            for i in range(0, len(self.adn), tam_okazaki)
        ]

        self.fragmentos_okazaki = []

        for molde in reversed(fragmentos_molde):
            nuevo = molde.reverse_complement()
            n = min(tam_cebador, len(nuevo))
            cebador = nuevo[:n].transcribe()

            self.fragmentos_okazaki.append(
                {
                    "fragmento": nuevo,
                    "cebador": cebador,
                }
            )

        self.cadena_rezagada = self.adn.reverse_complement()
        self.replicacion_realizada = True

    def mostrar_replicacion(self):
        if not self.comprobar_adn():
            return

        self.calcular_replicacion()
        self.marcar_activo(self.btn_rep)
        self.limpiar()
        cont = self.nuevo_contenedor()

        self.encabezado(
            cont,
            "Replicación del ADN",
            "Síntesis de dos moléculas hijas mediante replicación semiconservativa.",
        )
        self.paso_actual(cont, 2)

        metricas = ctk.CTkFrame(cont, fg_color="transparent")
        metricas.pack(fill="x", pady=(0, 16))
        for i in range(3):
            metricas.grid_columnconfigure(i, weight=1)

        self.metrica(metricas, 0, "Hebra líder", "Continua", "Síntesis 5' → 3'", SUCCESS)
        self.metrica(
            metricas,
            1,
            "Hebra rezagada",
            f"{len(self.fragmentos_okazaki)} fragmentos",
            "Fragmentos de Okazaki",
            WARNING,
        )
        self.metrica(metricas, 2, "Resultado", "2 moléculas", "Replicación semiconservativa", DNA)

        tabs = ctk.CTkTabview(
            cont,
            fg_color=PANEL,
            segmented_button_fg_color=PANEL_2,
            segmented_button_selected_color="#246DA8",
            segmented_button_selected_hover_color="#2F80C2",
            border_width=1,
            border_color=BORDER,
            corner_radius=13,
        )
        tabs.pack(fill="x", pady=(0, 16))

        for nombre in ["Esquema", "Secuencias", "Okazaki", "Enzimas"]:
            tabs.add(nombre)

        self._tab_replicacion_esquema(tabs.tab("Esquema"))
        self._tab_replicacion_secuencias(tabs.tab("Secuencias"))
        self._tab_replicacion_okazaki(tabs.tab("Okazaki"))
        self._tab_replicacion_enzimas(tabs.tab("Enzimas"))

        resultado = self.panel(cont, "Resultado de la replicación")

        hijas = ctk.CTkFrame(resultado, fg_color="transparent")
        hijas.pack(fill="x", padx=22, pady=(7, 18))
        hijas.grid_columnconfigure(0, weight=1)
        hijas.grid_columnconfigure(1, weight=1)

        self._molecula_hija(hijas, 0, "Molécula hija 1", "parental", "nueva")
        self._molecula_hija(hijas, 1, "Molécula hija 2", "nueva", "parental")

        ctk.CTkLabel(
            resultado,
            text="Cada molécula hija contiene una hebra parental y una hebra de nueva síntesis.",
            text_color=MUTED,
            font=ctk.CTkFont(size=12),
        ).pack(anchor="w", padx=22, pady=(0, 7))

        ctk.CTkLabel(
            resultado,
            text="Nota: se usan fragmentos de 500 nt como simplificación didáctica para visualizar la cadena rezagada.",
            text_color=SUBTLE,
            font=ctk.CTkFont(size=10),
        ).pack(anchor="w", padx=22, pady=(0, 18))

        self.boton_siguiente(cont, "Continuar a transcripción  →", self.mostrar_transcripcion)

    def _tab_replicacion_esquema(self, tab):
        canvas = tk.Canvas(
            tab,
            height=260,
            bg=PANEL,
            highlightthickness=0,
        )
        canvas.pack(fill="x", padx=18, pady=16)
        canvas.bind("<Configure>", lambda e: self._dibujar_rep_compacto(canvas, e.width))

    def _dibujar_rep_compacto(self, canvas, ancho):
        canvas.delete("all")
        ancho = max(ancho, 760)

        x0 = 40
        xf = ancho - 40
        xfork = int(ancho * 0.57)

        y1 = 95
        y2 = 155

        # El ADN sin separar queda a la derecha.
        # Por tanto, la horquilla avanza hacia la derecha.
        #
        # Hebra parental superior: 5' -> 3' = molde rezagado
        # Hebra parental inferior: 3' -> 5' = molde líder

        # ====================================================
        # HEBRAS PARENTALES ABIERTAS
        # ====================================================

        canvas.create_line(
            x0,
            45,
            xfork,
            y1,
            fill=DNA,
            width=4
        )

        canvas.create_line(
            x0,
            205,
            xfork,
            y2,
            fill=DNA,
            width=4
        )

        # ====================================================
        # ADN DOBLE TODAVÍA SIN SEPARAR
        # ====================================================

        canvas.create_line(
            xfork,
            y1,
            xf,
            y1,
            fill=DNA,
            width=4
        )

        canvas.create_line(
            xfork,
            y2,
            xf,
            y2,
            fill=DNA,
            width=4
        )

        for x in range(
            xfork + 14,
            xf - 8,
            24
        ):
            canvas.create_line(
                x,
                y1 + 6,
                x,
                y2 - 6,
                fill="#38404A",
                width=1
            )

        # ====================================================
        # AVANCE DE LA HORQUILLA
        # ====================================================

        canvas.create_text(
            xfork,
            18,
            text="avance de la horquilla",
            fill=MUTED,
            font=("Arial", 9)
        )

        canvas.create_line(
            xfork - 70,
            36,
            xfork + 70,
            36,
            fill=MUTED,
            width=2,
            arrow=tk.LAST
        )

        # ====================================================
        # CADENA LÍDER
        # ====================================================
        # Se sintetiza continuamente 5' -> 3'
        # en el mismo sentido que avanza la horquilla.

        canvas.create_line(
            x0 + 25,
            197,
            xfork - 16,
            y2 + 4,
            fill=SUCCESS,
            width=5,
            arrow=tk.LAST
        )

        # ====================================================
        # CADENA REZAGADA
        # ====================================================
        # Los fragmentos de Okazaki también se sintetizan
        # 5' -> 3', pero alejándose de la horquilla.

        span = xfork - x0

        def y_superior(x):
            fraccion = (x - x0) / span
            return 45 + (y1 - 45) * fraccion

        tramos = [
            (0.92, 0.68),
            (0.62, 0.38),
            (0.32, 0.08)
        ]

        for inicio_frac, fin_frac in tramos:

            x_inicio = x0 + span * inicio_frac
            x_fin = x0 + span * fin_frac

            y_inicio = y_superior(x_inicio) + 9
            y_fin = y_superior(x_fin) + 9

            canvas.create_line(
                x_inicio,
                y_inicio,
                x_fin,
                y_fin,
                fill=WARNING,
                width=5,
                arrow=tk.LAST
            )

        # ====================================================
        # HELICASA
        # ====================================================

        canvas.create_oval(
            xfork - 22,
            103,
            xfork + 22,
            147,
            fill="#8B72D8",
            outline=""
        )

        canvas.create_text(
            xfork,
            125,
            text="H",
            fill=TEXT,
            font=("Arial", 12, "bold")
        )

        canvas.create_text(
            xfork,
            168,
            text="helicasa",
            fill=MUTED,
            font=("Arial", 9)
        )

        # ====================================================
        # POLARIDADES
        # ====================================================

        # Hebra parental superior: 5' -> 3'

        canvas.create_text(
            x0,
            28,
            text="5'",
            fill=MUTED,
            anchor="w",
            font=("Arial", 10)
        )

        canvas.create_text(
            xf,
            78,
            text="3'",
            fill=MUTED,
            anchor="e",
            font=("Arial", 10)
        )

        # Hebra parental inferior: 3' -> 5'

        canvas.create_text(
            x0,
            222,
            text="3'",
            fill=MUTED,
            anchor="w",
            font=("Arial", 10)
        )

        canvas.create_text(
            xf,
            178,
            text="5'",
            fill=MUTED,
            anchor="e",
            font=("Arial", 10)
        )

        # ====================================================
        # ETIQUETAS
        # ====================================================

        canvas.create_text(
            int(ancho * 0.28),
            18,
            text="Cadena rezagada · síntesis discontinua",
            fill=WARNING,
            font=("Arial", 11, "bold")
        )

        canvas.create_text(
            int(ancho * 0.28),
            228,
            text="Cadena líder · síntesis continua 5' → 3'",
            fill=SUCCESS,
            font=("Arial", 11, "bold")
        )

        # ====================================================
        # LEYENDA
        # ====================================================

        canvas.create_text(
            int(ancho * 0.54),
            238,
            text=(
                "Azul: hebra parental     "
                "Verde: nueva hebra líder     "
                "Amarillo: fragmentos de Okazaki"
            ),
            fill=SUBTLE,
            anchor="w",
            font=("Arial", 9)
        )

    def _tab_replicacion_secuencias(self, tab):
        n = min(68, len(self.adn))
        molde_lider = self.adn[:n].complement()

        texto = (
            f"Hebra molde líder     3'  {molde_lider}  5'\n"
            f"Nueva hebra líder     5'  {self.cadena_lider[:n]}  3'\n\n"
            f"Cebador ARN           5'  {self.cebador_lider}  3'\n\n"
            f"Cadena rezagada final 5'  {self.cadena_rezagada[:n]}  3'"
        )

        caja = ctk.CTkTextbox(
            tab,
            height=210,
            font=(MONO, 12),
            fg_color="#14171A",
            text_color="#DDE2E7",
            border_width=1,
            border_color=BORDER,
            corner_radius=9,
            wrap="none",
        )
        caja.pack(fill="x", padx=18, pady=18)
        caja.insert("1.0", texto)
        caja.configure(state="disabled")

    def _tab_replicacion_okazaki(self, tab):
        lineas = []

        for i, datos in enumerate(self.fragmentos_okazaki, start=1):
            frag = datos["fragmento"]
            ceb = datos["cebador"]
            lineas.append(
                f"{i:>2}. {len(frag):>4} nt   cebador 5' {ceb} 3'   "
                f"ADN 5' {self.formatear(frag, 42)} 3'"
            )

        texto = (
            "Fragmentos simulados de la cadena rezagada\n"
            "──────────────────────────────────────────\n"
            + "\n".join(lineas)
        )

        caja = ctk.CTkTextbox(
            tab,
            height=250,
            font=(MONO, 11),
            fg_color="#14171A",
            text_color="#DDE2E7",
            border_width=1,
            border_color=BORDER,
            corner_radius=9,
            wrap="none",
        )
        caja.pack(fill="x", padx=18, pady=18)
        caja.insert("1.0", texto)
        caja.configure(state="disabled")

    def _tab_replicacion_enzimas(self, tab):
        elementos = [
            ("Helicasa", "Separa las dos hebras de ADN."),
            ("Proteínas SSB", "Estabilizan las hebras abiertas."),
            ("Primasa", "Sintetiza cebadores de ARN."),
            ("ADN polimerasa III", "Elonga las nuevas hebras 5' → 3' en E. coli."),
            ("ADN polimerasa I", "Reemplaza los cebadores de ARN por ADN en E. coli."),
            ("Ligasa", "Sella los fragmentos de Okazaki."),
        ]

        grid = ctk.CTkFrame(tab, fg_color="transparent")
        grid.pack(fill="x", padx=14, pady=14)

        for col in range(2):
            grid.grid_columnconfigure(col, weight=1)

        for i, (nombre, desc) in enumerate(elementos):
            card = ctk.CTkFrame(
                grid,
                fg_color=PANEL_2,
                corner_radius=9,
                border_width=1,
                border_color=BORDER,
            )
            card.grid(row=i // 2, column=i % 2, sticky="nsew", padx=5, pady=5)

            ctk.CTkLabel(
                card,
                text=nombre,
                text_color=TEXT,
                font=ctk.CTkFont(size=13, weight="bold"),
            ).pack(anchor="w", padx=14, pady=(11, 2))

            ctk.CTkLabel(
                card,
                text=desc,
                text_color=MUTED,
                justify="left",
                wraplength=380,
                font=ctk.CTkFont(size=11),
            ).pack(anchor="w", padx=14, pady=(2, 11))

    def _molecula_hija(self, padre, columna, titulo, superior_tipo, inferior_tipo):
        frame = ctk.CTkFrame(
            padre,
            fg_color=PANEL_2,
            corner_radius=9,
            border_width=1,
            border_color=BORDER,
        )
        frame.grid(row=0, column=columna, sticky="nsew", padx=5)

        ctk.CTkLabel(
            frame,
            text=titulo,
            text_color=TEXT,
            font=ctk.CTkFont(size=13, weight="bold"),
        ).pack(anchor="w", padx=14, pady=(12, 7))

        color_sup = DNA if superior_tipo == "parental" else SUCCESS
        color_inf = DNA if inferior_tipo == "parental" else SUCCESS

        ctk.CTkLabel(
            frame,
            text=f"● Hebra {superior_tipo}",
            text_color=color_sup,
            font=ctk.CTkFont(size=11, weight="bold"),
        ).pack(anchor="w", padx=14)

        ctk.CTkLabel(
            frame,
            text=f"5' {self.formatear(self.adn, 42)} 3'",
            text_color=MUTED,
            font=(MONO, 9),
        ).pack(anchor="w", padx=14, pady=(2, 7))

        ctk.CTkLabel(
            frame,
            text=f"● Hebra {inferior_tipo}",
            text_color=color_inf,
            font=ctk.CTkFont(size=11, weight="bold"),
        ).pack(anchor="w", padx=14)

        ctk.CTkLabel(
            frame,
            text=f"3' {self.formatear(self.adn.complement(), 42)} 5'",
            text_color=MUTED,
            font=(MONO, 9),
        ).pack(anchor="w", padx=14, pady=(2, 12))

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
        self.marcar_activo(self.btn_trans)
        self.limpiar()
        cont = self.nuevo_contenedor()

        self.encabezado(
            cont,
            "Transcripción",
            "La ARN polimerasa utiliza la hebra molde para sintetizar ARN mensajero.",
        )
        self.paso_actual(cont, 3)

        metricas = ctk.CTkFrame(cont, fg_color="transparent")
        metricas.pack(fill="x", pady=(0, 16))
        for i in range(3):
            metricas.grid_columnconfigure(i, weight=1)

        self.metrica(metricas, 0, "Entrada", f"{len(self.adn)} nt", "ADN codificante", DNA)
        self.metrica(metricas, 1, "Enzima", "ARN polimerasa", "Lee la hebra molde", TEXT)
        self.metrica(metricas, 2, "Salida", f"{len(self.arn_m)} nt", "ARN mensajero", RNA)

        tabs = ctk.CTkTabview(
            cont,
            fg_color=PANEL,
            segmented_button_fg_color=PANEL_2,
            segmented_button_selected_color="#246DA8",
            segmented_button_selected_hover_color="#2F80C2",
            border_width=1,
            border_color=BORDER,
            corner_radius=13,
        )
        tabs.pack(fill="x", pady=(0, 16))

        for nombre in ["Esquema", "Secuencias", "Reglas"]:
            tabs.add(nombre)

        self._tab_trans_esquema(tabs.tab("Esquema"))
        self._tab_trans_secuencias(tabs.tab("Secuencias"))
        self._tab_trans_reglas(tabs.tab("Reglas"))

        self.boton_siguiente(cont, "Continuar a traducción  →", self.mostrar_traduccion)

    def _tab_trans_esquema(self, tab):
        canvas = tk.Canvas(
            tab,
            height=230,
            bg=PANEL,
            highlightthickness=0,
        )
        canvas.pack(fill="x", padx=18, pady=16)
        canvas.bind("<Configure>", lambda e: self._dibujar_trans_compacto(canvas, e.width))

    def _dibujar_trans_compacto(self, canvas, ancho):
        canvas.delete("all")
        ancho = max(ancho, 760)

        y = 105
        x1 = int(ancho * 0.18)
        x2 = int(ancho * 0.50)
        x3 = int(ancho * 0.82)

        # Cajas
        self._canvas_box(canvas, x1, y, 180, 70, "ADN molde", "3' → 5'", DNA)
        self._canvas_box(canvas, x2, y, 210, 70, "ARN polimerasa", "síntesis", "#8B72D8")
        self._canvas_box(canvas, x3, y, 180, 70, "ARNm", "5' → 3'", RNA)

        canvas.create_line(x1 + 105, y, x2 - 120, y, fill="#56606A", width=2, arrow=tk.LAST)
        canvas.create_line(x2 + 120, y, x3 - 105, y, fill="#56606A", width=2, arrow=tk.LAST)

        canvas.create_text(
            ancho // 2,
            188,
            text="La secuencia del ARNm coincide con la hebra codificante salvo que U sustituye a T.",
            fill=MUTED,
            font=("Arial", 10),
        )

    def _canvas_box(self, canvas, x, y, w, h, titulo, sub, color):
        canvas.create_rectangle(
            x - w // 2,
            y - h // 2,
            x + w // 2,
            y + h // 2,
            fill="#20252A",
            outline=color,
            width=2,
        )
        canvas.create_text(
            x,
            y - 8,
            text=titulo,
            fill=TEXT,
            font=("Arial", 12, "bold"),
        )
        canvas.create_text(
            x,
            y + 14,
            text=sub,
            fill=color,
            font=("Arial", 9),
        )

    def _tab_trans_secuencias(self, tab):
        n = min(72, len(self.adn))

        texto = (
            f"Codificante   5'  {self.adn[:n]}  3'\n"
            f"Molde         3'  {self.adn[:n].complement()}  5'\n"
            f"                  {'│' * n}\n"
            f"ARNm          5'  {self.arn_m[:n]}  3'\n\n"
            f"Longitud del ARNm: {len(self.arn_m)} nucleótidos"
        )

        caja = ctk.CTkTextbox(
            tab,
            height=210,
            font=(MONO, 12),
            fg_color="#14171A",
            text_color="#DDE2E7",
            border_width=1,
            border_color=BORDER,
            corner_radius=9,
            wrap="none",
        )
        caja.pack(fill="x", padx=18, pady=18)
        caja.insert("1.0", texto)
        caja.configure(state="disabled")

    def _tab_trans_reglas(self, tab):
        cards = ctk.CTkFrame(tab, fg_color="transparent")
        cards.pack(fill="x", padx=14, pady=14)
        cards.grid_columnconfigure(0, weight=1)
        cards.grid_columnconfigure(1, weight=1)

        datos = [
            ("Dirección de lectura", "La ARN polimerasa lee la hebra molde 3' → 5'."),
            ("Dirección de síntesis", "El ARNm se sintetiza 5' → 3'."),
            ("Complementariedad", "A ↔ U y C ↔ G durante la formación del ARN."),
            ("Diferencia ADN / ARN", "El ARN contiene uracilo (U) en lugar de timina (T)."),
        ]

        for i, (t, d) in enumerate(datos):
            card = ctk.CTkFrame(
                cards,
                fg_color=PANEL_2,
                corner_radius=9,
                border_width=1,
                border_color=BORDER,
            )
            card.grid(row=i // 2, column=i % 2, sticky="nsew", padx=5, pady=5)

            ctk.CTkLabel(
                card,
                text=t,
                text_color=TEXT,
                font=ctk.CTkFont(size=13, weight="bold"),
            ).pack(anchor="w", padx=14, pady=(11, 2))

            ctk.CTkLabel(
                card,
                text=d,
                text_color=MUTED,
                justify="left",
                wraplength=400,
                font=ctk.CTkFont(size=11),
            ).pack(anchor="w", padx=14, pady=(2, 11))

    # ========================================================
    # TRADUCCIÓN
    # ========================================================

    def calcular_traduccion(self):
        if self.traduccion_realizada:
            return

        self.calcular_transcripcion()

        rna = str(self.arn_m)
        inicio = rna.find("AUG")

        if inicio == -1:
            self.inicio_traduccion = None
            self.proteina = Seq("")
            self.traduccion_realizada = True
            return

        self.inicio_traduccion = inicio

        stops = {"UAA", "UAG", "UGA"}
        stop_pos = None
        stop_codon = None

        for pos in range(inicio, len(rna) - 2, 3):
            codon = rna[pos:pos + 3]

            if codon in stops:
                stop_pos = pos
                stop_codon = codon
                break

        if stop_pos is not None:
            region = self.arn_m[inicio:stop_pos]
            self.codon_stop = stop_codon
            self.posicion_stop_nt = stop_pos
            self.posicion_stop_codon = ((stop_pos - inicio) // 3) + 1
        else:
            region = self.arn_m[inicio:]
            resto = len(region) % 3

            if resto:
                region = region[:-resto]

            self.codon_stop = None
            self.posicion_stop_nt = None
            self.posicion_stop_codon = None

        self.proteina = region.translate(table=11) if len(region) >= 3 else Seq("")

        self.codones_traducidos = []

        for pos in range(0, len(region), 3):
            codon = region[pos:pos + 3]

            if len(codon) < 3:
                continue

            aa = codon.translate(table=11)

            # IMPORTANTE: anticodón de ARN, no ADN.
            # Se presenta antiparalelo al codón: 3' UAC 5' frente a 5' AUG 3'.
            anticodon = codon.complement_rna()

            self.codones_traducidos.append(
                {
                    "codon": str(codon),
                    "anticodon": str(anticodon),
                    "aa": str(aa),
                }
            )

        self.traduccion_realizada = True

    def mostrar_traduccion(self):
        if not self.comprobar_adn():
            return

        self.calcular_traduccion()
        self.marcar_activo(self.btn_trad)
        self.limpiar()
        cont = self.nuevo_contenedor()

        self.encabezado(
            cont,
            "Traducción",
            "El ribosoma lee los codones del ARNm y construye una secuencia de aminoácidos.",
        )
        self.paso_actual(cont, 4)

        if self.inicio_traduccion is None:
            aviso = self.panel(
                cont,
                "No se ha encontrado un codón AUG",
                "El simulador utiliza AUG como codón de inicio didáctico. Sin él no se inicia la traducción.",
            )
            ctk.CTkLabel(
                aviso,
                text="Prueba con otra secuencia o utiliza el ejemplo incluido.",
                text_color=WARNING,
                font=ctk.CTkFont(size=13, weight="bold"),
            ).pack(anchor="w", padx=22, pady=(4, 18))
            return

        metricas = ctk.CTkFrame(cont, fg_color="transparent")
        metricas.pack(fill="x", pady=(0, 16))
        for i in range(3):
            metricas.grid_columnconfigure(i, weight=1)

        self.metrica(
            metricas,
            0,
            "Inicio",
            f"AUG · nt {self.inicio_traduccion + 1}",
            "Codón de inicio",
            RNA,
        )
        self.metrica(
            metricas,
            1,
            "Proteína",
            f"{len(self.proteina)} aa",
            "Tabla genética bacteriana 11",
            PROTEIN,
        )
        self.metrica(
            metricas,
            2,
            "Terminación",
            self.codon_stop or "Sin STOP",
            "UAA · UAG · UGA",
            SUCCESS if self.codon_stop else WARNING,
        )

        tabs = ctk.CTkTabview(
            cont,
            fg_color=PANEL,
            segmented_button_fg_color=PANEL_2,
            segmented_button_selected_color="#246DA8",
            segmented_button_selected_hover_color="#2F80C2",
            border_width=1,
            border_color=BORDER,
            corner_radius=13,
        )
        tabs.pack(fill="x", pady=(0, 16))

        for nombre in ["Esquema", "Codones", "ARNt", "Proteína"]:
            tabs.add(nombre)

        self._tab_trad_esquema(tabs.tab("Esquema"))
        self._tab_trad_codones(tabs.tab("Codones"))
        self._tab_trad_arnt(tabs.tab("ARNt"))
        self._tab_trad_proteina(tabs.tab("Proteína"))

        nota = self.panel(
            cont,
            "Criterio utilizado",
            "Para secuencias genéricas se busca el primer AUG y se continúa en ese marco hasta UAA, UAG o UGA. Para lacZ de E. coli se utiliza la tabla genética bacteriana 11.",
        )

        ctk.CTkLabel(
            nota,
            text=(
                f"Resultado: {len(self.proteina)} aminoácidos"
                + (f" · STOP {self.codon_stop}" if self.codon_stop else " · no se encontró STOP en el marco")
            ),
            text_color=PROTEIN,
            font=ctk.CTkFont(size=13, weight="bold"),
        ).pack(anchor="w", padx=22, pady=(4, 18))

        self.boton_siguiente(cont, "Ver resultados finales  →", self.mostrar_resultados)

    def _tab_trad_esquema(self, tab):
        canvas = tk.Canvas(
            tab,
            height=245,
            bg=PANEL,
            highlightthickness=0,
        )
        canvas.pack(fill="x", padx=18, pady=16)
        canvas.bind("<Configure>", lambda e: self._dibujar_trad_compacto(canvas, e.width))

    def _dibujar_trad_compacto(self, canvas, ancho):
        canvas.delete("all")
        ancho = max(ancho, 760)

        datos = self.codones_traducidos[:5]
        if not datos:
            return

        centro = ancho // 2
        y_rna = 145

        # Línea ARNm
        canvas.create_line(
            80,
            y_rna,
            ancho - 80,
            y_rna,
            fill=RNA,
            width=4,
            arrow=tk.LAST,
        )

        canvas.create_text(
            80,
            y_rna + 24,
            text="5' ARNm",
            fill=RNA,
            anchor="w",
            font=("Arial", 10, "bold"),
        )

        # Ribosoma, limpio y pequeño
        canvas.create_rectangle(
            centro - 250,
            55,
            centro + 250,
            120,
            fill="#2A2536",
            outline="#8B72D8",
            width=2,
        )
        canvas.create_text(
            centro,
            75,
            text="RIBOSOMA",
            fill=TEXT,
            font=("Arial", 11, "bold"),
        )
        canvas.create_text(
            centro,
            98,
            text="lectura del ARNm en codones de tres nucleótidos",
            fill=MUTED,
            font=("Arial", 9),
        )

        separacion = 108
        inicio_x = centro - ((len(datos) - 1) * separacion) / 2

        for i, dato in enumerate(datos):
            x = inicio_x + i * separacion
            codon = dato["codon"]
            aa = dato["aa"]

            canvas.create_rectangle(
                x - 35,
                y_rna - 16,
                x + 35,
                y_rna + 16,
                fill="#173027",
                outline=RNA,
                width=1,
            )
            canvas.create_text(
                x,
                y_rna,
                text=codon,
                fill=TEXT,
                font=(MONO, 10, "bold"),
            )

            canvas.create_line(
                x,
                y_rna - 18,
                x,
                121,
                fill="#56606A",
                width=1,
            )

            canvas.create_text(
                x,
                188,
                text=aa,
                fill=PROTEIN,
                font=("Arial", 15, "bold"),
            )

        canvas.create_text(
            centro,
            218,
            text="Los aminoácidos se enlazan para formar la cadena polipeptídica.",
            fill=MUTED,
            font=("Arial", 10),
        )

    def _tab_trad_codones(self, tab):
        grid = ctk.CTkFrame(tab, fg_color="transparent")
        grid.pack(fill="x", padx=14, pady=14)

        for col in range(6):
            grid.grid_columnconfigure(col, weight=1)

        for i, dato in enumerate(self.codones_traducidos[:12]):
            card = ctk.CTkFrame(
                grid,
                fg_color=PANEL_2,
                corner_radius=8,
                border_width=1,
                border_color=BORDER,
            )
            card.grid(row=i // 6, column=i % 6, sticky="nsew", padx=4, pady=4)

            aa = dato["aa"]
            nombre = AMINOACIDOS.get(aa, aa)

            ctk.CTkLabel(
                card,
                text=dato["codon"],
                text_color=RNA,
                font=ctk.CTkFont(size=14, weight="bold"),
            ).pack(pady=(10, 2))

            ctk.CTkLabel(
                card,
                text=aa,
                text_color=PROTEIN,
                font=ctk.CTkFont(size=18, weight="bold"),
            ).pack()

            ctk.CTkLabel(
                card,
                text=nombre,
                text_color=SUBTLE,
                font=ctk.CTkFont(size=9),
            ).pack(pady=(2, 10))

    def _tab_trad_arnt(self, tab):
        lineas = [
            "El anticodón se muestra antiparalelo al codón del ARNm.",
            "",
            "ARNm                 ARNt",
            "────────────────────────────────────",
        ]

        for i, dato in enumerate(self.codones_traducidos[:10], start=1):
            lineas.append(
                f"{i:>2}.  5' {dato['codon']} 3'   ↔   3' {dato['anticodon']} 5'   "
                f"→ {dato['aa']} ({AMINOACIDOS.get(dato['aa'], dato['aa'])})"
            )

        texto = "\n".join(lineas)

        caja = ctk.CTkTextbox(
            tab,
            height=250,
            font=(MONO, 11),
            fg_color="#14171A",
            text_color="#DDE2E7",
            border_width=1,
            border_color=BORDER,
            corner_radius=9,
            wrap="none",
        )
        caja.pack(fill="x", padx=18, pady=18)
        caja.insert("1.0", texto)
        caja.configure(state="disabled")

    def _tab_trad_proteina(self, tab):
        if len(self.proteina) == 0:
            texto = "No se ha sintetizado proteína."
        else:
            texto = (
                f"Longitud: {len(self.proteina)} aminoácidos\n\n"
                f"Primeros 60:\n{self.proteina[:60]}\n\n"
                f"Últimos 60:\n{self.proteina[-60:]}\n\n"
                f"Secuencia completa:\n{self.proteina}"
            )

        caja = ctk.CTkTextbox(
            tab,
            height=270,
            font=(MONO, 11),
            fg_color="#14171A",
            text_color="#DDE2E7",
            border_width=1,
            border_color=BORDER,
            corner_radius=9,
            wrap="word",
        )
        caja.pack(fill="x", padx=18, pady=18)
        caja.insert("1.0", texto)
        caja.configure(state="disabled")

    # ========================================================
    # RESULTADOS
    # ========================================================

    def mostrar_resultados(self):
        if not self.comprobar_adn():
            return

        self.calcular_replicacion()
        self.calcular_transcripcion()
        self.calcular_traduccion()

        self.marcar_activo(self.btn_res)
        self.limpiar()
        cont = self.nuevo_contenedor()

        self.encabezado(
            cont,
            "Resultado de la simulación",
            "Resumen compacto del flujo completo de información genética.",
        )

        metricas = ctk.CTkFrame(cont, fg_color="transparent")
        metricas.pack(fill="x", pady=(0, 16))

        for i in range(3):
            metricas.grid_columnconfigure(i, weight=1)

        self.metrica(metricas, 0, "ADN", f"{len(self.adn)} nt", "Secuencia inicial", DNA)
        self.metrica(metricas, 1, "ARNm", f"{len(self.arn_m)} nt", "Producto de transcripción", RNA)

        prot = (
            f"{len(self.proteina)} aa"
            if self.inicio_traduccion is not None
            else "No traducida"
        )
        self.metrica(metricas, 2, "Proteína", prot, "Producto de traducción", PROTEIN)

        flujo = self.panel(cont, "Dogma central")

        barra = ctk.CTkFrame(flujo, fg_color="transparent")
        barra.pack(fill="x", padx=22, pady=(8, 18))

        pasos = [
            ("ADN", DNA),
            ("Replicación", TEXT),
            ("ADN", DNA),
            ("Transcripción", TEXT),
            ("ARNm", RNA),
            ("Traducción", TEXT),
            ("Proteína", PROTEIN),
        ]

        for i, (nombre, color) in enumerate(pasos):
            ctk.CTkLabel(
                barra,
                text=nombre,
                text_color=color,
                font=ctk.CTkFont(size=14, weight="bold"),
            ).pack(side="left")

            if i < len(pasos) - 1:
                ctk.CTkLabel(
                    barra,
                    text="  →  ",
                    text_color=SUBTLE,
                    font=ctk.CTkFont(size=14),
                ).pack(side="left")

        resumen = self.panel(cont, "Resumen")

        inicio_txt = (
            str(self.inicio_traduccion + 1)
            if self.inicio_traduccion is not None
            else "No encontrado"
        )

        stop_txt = self.codon_stop or "No encontrado"

        filas = [
            ("Secuencia", self.registro.id),
            ("Longitud del ADN", f"{len(self.adn)} nt"),
            ("Fragmentos de Okazaki simulados", str(len(self.fragmentos_okazaki))),
            ("Longitud del ARNm", f"{len(self.arn_m)} nt"),
            ("Primer AUG", f"posición {inicio_txt}"),
            ("Proteína", prot),
            ("Codón STOP", stop_txt),
        ]

        tabla = ctk.CTkFrame(resumen, fg_color="transparent")
        tabla.pack(fill="x", padx=22, pady=(7, 17))
        tabla.grid_columnconfigure(1, weight=1)

        for fila, (campo, valor) in enumerate(filas):
            ctk.CTkLabel(
                tabla,
                text=campo,
                text_color=MUTED,
                font=ctk.CTkFont(size=12),
            ).grid(row=fila, column=0, sticky="w", padx=(0, 28), pady=5)

            ctk.CTkLabel(
                tabla,
                text=valor,
                text_color=TEXT,
                font=ctk.CTkFont(size=12, weight="bold"),
            ).grid(row=fila, column=1, sticky="w", pady=5)

        final = self.panel(
            cont,
            "Simulación completada",
            "El programa ha representado ADN → ADN, ADN → ARN y ARN → proteína mediante operaciones de Biopython.",
        )

        ctk.CTkLabel(
            final,
            text="✓ Flujo completo procesado",
            text_color=SUCCESS,
            font=ctk.CTkFont(size=14, weight="bold"),
        ).pack(anchor="w", padx=22, pady=(3, 18))

        acciones = ctk.CTkFrame(cont, fg_color="transparent")
        acciones.pack(fill="x", pady=(0, 18))

        ctk.CTkButton(
            acciones,
            text="Volver al inicio",
            command=self.mostrar_inicio,
            height=42,
            width=165,
            corner_radius=8,
        ).pack(side="left", padx=(0, 8))

        ctk.CTkButton(
            acciones,
            text="Cargar otra secuencia",
            command=self.cargar_archivo,
            height=42,
            width=185,
            corner_radius=8,
            fg_color=PANEL_2,
            hover_color="#30363D",
        ).pack(side="left")

    # ========================================================
    # UTILIDADES
    # ========================================================

    def comprobar_adn(self):
        if self.adn is None:
            self.mostrar_error("Primero debes cargar una secuencia de ADN.")
            return False
        return True

    def obtener_nombre_archivo(self):
        if self.ruta_archivo is not None:
            return Path(self.ruta_archivo).name

        if self.registro is not None:
            return f"Secuencia activa: {self.registro.id}"

        return "Ninguna secuencia cargada"

    def mostrar_error(self, mensaje):
        ventana = ctk.CTkToplevel(self)
        ventana.title("Aviso")
        ventana.geometry("440x220")
        ventana.resizable(False, False)
        ventana.transient(self)
        ventana.grab_set()

        ctk.CTkLabel(
            ventana,
            text="Aviso",
            text_color=TEXT,
            font=ctk.CTkFont(size=21, weight="bold"),
        ).pack(pady=(28, 6))

        ctk.CTkLabel(
            ventana,
            text=mensaje,
            text_color=MUTED,
            wraplength=360,
            justify="center",
        ).pack(padx=25, pady=10)

        ctk.CTkButton(
            ventana,
            text="Aceptar",
            command=ventana.destroy,
            width=120,
        ).pack(pady=10)


if __name__ == "__main__":
    app = SimuladorDogmaCentral()
    app.mainloop()
