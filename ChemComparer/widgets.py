# біблиотеки
from PyQt5.QtWidgets import (QWidget, QFileDialog, QVBoxLayout, QHBoxLayout, QRadioButton, QLineEdit, QPushButton, QLabel, QTableWidget, QTableWidgetItem, QSizePolicy, QFrame, QMenuBar, QAction, QMessageBox, QMenu, QDialog, QScrollArea, QGridLayout, QStatusBar, QTabBar)
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPixmap, QIcon
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors, rdMolDescriptors
import py3Dmol
from io import BytesIO
import json
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import pandas as pd

class ComparisonDiagramTab(QWidget):
    def __init__(self, molecules_data, parent=None):
        super().__init__(parent)
        self.molecules_data = molecules_data
        self.parent = parent
        
        layout = QVBoxLayout()
        
        # кнопка запуску діаграми
        button_layout = QHBoxLayout()
        lang = self.parent.current_language
        window_button_text = "Open diagram in window" if lang == 'en' else "Diagramm in Fenster öffnen" if lang == 'de' else "Відкрити діаграму у вікні" if lang == 'uk' else "Ouvrir le diagramme dans la fenêtre" if lang == 'fr' else "Open diagram in window"
        self.window_button = QPushButton(window_button_text)
        self.window_button.clicked.connect(self.open_diagram_window)
        button_layout.addWidget(self.window_button)
        button_layout.addStretch()
        layout.addLayout(button_layout)
        
        self.figure = Figure(figsize=(10, 8), dpi=100)
        self.canvas = FigureCanvas(self.figure)
        layout.addWidget(self.canvas)
        

        self.create_comparison_diagram()
        
        self.setLayout(layout)
    
    def create_comparison_diagram(self):
        self.figure.clear()
        
        # створення діаграми
        lang = self.parent.current_language
        if lang == 'en':
            properties = [
                'Molecular weight', 'LogP', 'Number of atoms', 'Number of bonds',
                'Aromatic atoms', 'Rings', 'H donors', 'H acceptors'
            ]
        elif lang == 'de':
            properties = [
                'Molekulargewicht', 'LogP', 'Anzahl der Atome', 'Anzahl der Bindungen',
                'Aromatische Atome', 'Ringe', 'H-Spender', 'H-Akzeptoren'
            ]
        elif lang == 'uk':
            properties = [
                'Молекулярна вага', 'LogP', 'Число атомів', 'Число зв\'язків',
                'Ароматичні атоми', 'Цикли', 'Донори H', 'Акцептори H'
            ]
        elif lang == 'fr':
            properties = [
                'Poids moléculaire', 'LogP', 'Nombre d\'atomes', 'Nombre de liaisons',
                'Atomes aromatiques', 'Cycles', 'Donneurs H', 'Accepteurs H'
            ]
        
        # перевірка на мову
        if lang == 'en':
            mol_names = [f"Molecule {i+1}" for i in range(len(self.molecules_data))]
        elif lang == 'de':
            mol_names = [f"Molekül {i+1}" for i in range(len(self.molecules_data))]
        elif lang == 'uk':
            mol_names = [f"Молекула {i+1}" for i in range(len(self.molecules_data))]
        elif lang == 'fr':
            mol_names = [f"Molécule {i+1}" for i in range(len(self.molecules_data))]
        data = {}
        
        for prop in properties:
            data[prop] = [self.molecules_data[i].get(prop, 0) for i in range(len(self.molecules_data))]
        

        df = pd.DataFrame(data, index=mol_names)
        

        ax = self.figure.add_subplot(111)
        

        df.plot(kind='bar', ax=ax, rot=0)
        if lang == 'en':
            ax.set_title('Comparison of molecule properties')
            ax.set_ylabel('Values')
        elif lang == 'de':
            ax.set_title('Vergleich der Moleküleigenschaften')
            ax.set_ylabel('Werte')
        elif lang == 'uk':
            ax.set_title('Порівняння властивостей молекул')
            ax.set_ylabel('Значення')
        elif lang == 'fr':
            ax.set_title('Comparaison des propriétés des molécules')
            ax.set_ylabel('Valeurs')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        

        self.figure.tight_layout()
        

        self.canvas.draw()
    # пока не використовується
    def open_diagram_window(self):
        try:

            diagram_window = QMainWindow()
            lang = self.parent.current_language
            if lang == 'en':
                diagram_window.setWindowTitle("Comparative Diagram")
            elif lang == 'de':
                diagram_window.setWindowTitle("Vergleichsdiagramm")
            elif lang == 'uk':
                diagram_window.setWindowTitle("Порівняльна діаграма")
            elif lang == 'fr':
                diagram_window.setWindowTitle("Diagramme comparatif")
            diagram_window.setGeometry(100, 100, 1000, 700)
            
            figure = Figure(figsize=(10, 8), dpi=100)
            canvas = FigureCanvas(figure)
            
            lang = self.parent.current_language
            if lang == 'en':
                properties = [
                    'Molecular weight', 'LogP', 'Number of atoms', 'Number of bonds',
                    'Aromatic atoms', 'Rings', 'H donors', 'H acceptors'
                ]
            elif lang == 'de':
                properties = [
                    'Molekulargewicht', 'LogP', 'Anzahl der Atome', 'Anzahl der Bindungen',
                    'Aromatische Atome', 'Ringe', 'H-Spender', 'H-Akzeptoren'
                ]
            elif lang == 'uk':
                properties = [
                    'Молекулярна вага', 'LogP', 'Число атомів', 'Число зв\'язків',
                    'Ароматичні атоми', 'Цикли', 'Донори H', 'Акцептори H'
                ]
            elif lang == 'fr':
                properties = [
                    'Poids moléculaire', 'LogP', 'Nombre d\'atomes', 'Nombre de liaisons',
                    'Atomes aromatiques', 'Cycles', 'Donneurs H', 'Accepteurs H'
                ]
            
            ax = figure.add_subplot(111)
            if lang == 'en':
                mol_names = [f"Molecule {i+1}" for i in range(len(self.molecules_data))]
            elif lang == 'de':
                mol_names = [f"Molekül {i+1}" for i in range(len(self.molecules_data))]
            elif lang == 'uk':
                mol_names = [f"Молекула {i+1}" for i in range(len(self.molecules_data))]
            elif lang == 'fr':
                mol_names = [f"Molécule {i+1}" for i in range(len(self.molecules_data))]
            data = {}
            
            for prop in properties:
                data[prop] = [self.molecules_data[i].get(prop, 0) for i in range(len(self.molecules_data))]
            
            df = pd.DataFrame(data, index=mol_names)
            df.plot(kind='bar', ax=ax, rot=0)
            if lang == 'en':
                ax.set_title('Comparison of molecule properties')
                ax.set_ylabel('Values')
            elif lang == 'de':
                ax.set_title('Vergleich der Moleküleigenschaften')
                ax.set_ylabel('Werte')
            elif lang == 'uk':
                ax.set_title('Порівняння властивостей молекул')
                ax.set_ylabel('Значення')
            elif lang == 'fr':
                ax.set_title('Comparaison des propriétés des molécules')
                ax.set_ylabel('Valeurs')
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            figure.tight_layout()
            
            diagram_window.setCentralWidget(canvas)
            diagram_window.show()
        except Exception as e:
            lang = self.parent.current_language
            error_title = "Error" if lang == 'en' else "Fehler" if lang == 'de' else "Помилка" if lang == 'uk' else "Erreur" if lang == 'fr' else "Error"
            QMessageBox.warning(self, error_title, f"Error opening diagram in window: {str(e)}")

class MainTabWidget(QWidget):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent
        self.current_mol = None
        self.current_mol_2 = None
        self.current_mol_3 = None
        self.web_view = None
        self.current_mol_block = ""
        self.visualization_frame = None
        self.visualization_style = 'stick'
        self.show_hydrogens = True

        self.lang = self.parent.current_language

        self.main_tab_layout = QVBoxLayout(self)

        self.info_label = QLabel()
        self.info_label.setAlignment(Qt.AlignCenter)
        self.info_label.setStyleSheet("""
            QLabel {
                background-color: #f8f9fa;
                border: 2px dashed #dee2e6;
                border-radius: 10px;
                padding: 20px;
                margin: 20px;
                font-size: 14px;
                color: #6c757d;
            }
        """)
        if self.lang == 'en':
            info_text = "<h3>Welcome to ChemComparer 0.8.1!</h3><p>To start work, enter the molecule SMILES code in the field above.</p><p><strong>SMILES</strong> - is a text format for representing chemical structures.</p><p><em>Example: CCO for ethanol, c1ccccc1 for benzene</em></p><p>Use the 'Help' menu for additional information.</p>"
        elif self.lang == 'de':
            info_text = "<h3>Willkommen bei ChemComparer 0.8.1!</h3><p>Um zu starten, geben Sie den SMILES-Code des Moleküls in das Feld oben ein.</p><p><strong>SMILES</strong> - ist ein Textformat zur Darstellung chemischer Strukturen.</p><p><em>Beispiel: CCO für Ethanol, c1ccccc1 für Benzol</em></p><p>Verwenden Sie das Menü 'Hilfe' für zusätzliche Informationen.</p>"
        elif self.lang == 'uk':
            info_text = "<h3>Ласкаво просимо до ChemComparer 0.8.1!</h3><p>Для початку роботи введіть SMILES-код молекули в поле вище.</p><p><strong>SMILES</strong> - це текстовий формат для представлення хімічних структур.</p><p><em>Приклад: CCO для етанолу, c1ccccc1 для бензолу</em></p><p>Використовуйте меню 'Допомога' для додаткової інформації.</p>"
        elif self.lang == 'fr':
            info_text = "<h3>Bienvenue dans ChemComparer 0.8.1 !</h3><p>Pour commencer, entrez le code SMILES de la molécule dans le champ ci-dessus.</p><p><strong>SMILES</strong> - est un format texte pour représenter les structures chimiques.</p><p><em>Exemple : CCO pour l'éthanol, c1ccccc1 pour le benzène</em></p><p>Utilisez le menu 'Aide' pour des informations supplémentaires.</p>"
        self.info_label.setText(info_text)
        self.info_label.setWordWrap(True)


        input_layout = QHBoxLayout()
        self.line_edit = QLineEdit()
        if self.lang == 'en':
            placeholder_text = "Enter SMILES..."
        elif self.lang == 'de':
            placeholder_text = "SMILES eingeben..."
        elif self.lang == 'uk':
            placeholder_text = "Введіть SMILES..."
        elif self.lang == 'fr':
            placeholder_text = "Entrez SMILES..."
        self.line_edit.setPlaceholderText(placeholder_text)
        if self.lang == 'en':
            button_text = "Check SMILES and launch"
        elif self.lang == 'de':
            button_text = "SMILES überprüfen und starten"
        elif self.lang == 'uk':
            button_text = "Перевірити SMILES і запустити"
        elif self.lang == 'fr':
            button_text = "Vérifier SMILES et lancer"
        self.button = QPushButton(button_text)
        self.button.clicked.connect(self.check_smiles)
        
        if self.lang == 'en':
            radio_3d_text = "3D model"
        elif self.lang == 'de':
            radio_3d_text = "3D-Modell"
        elif self.lang == 'uk':
            radio_3d_text = "3D модель"
        elif self.lang == 'fr':
            radio_3d_text = "Modèle 3D"
        self.radio_3d = QRadioButton(radio_3d_text)
        if self.lang == 'en':
            radio_2d_text = "2D drawing"
        elif self.lang == 'de':
            radio_2d_text = "2D-Zeichnung"
        elif self.lang == 'uk':
            radio_2d_text = "2D малюнок"
        elif self.lang == 'fr':
            radio_2d_text = "Dessin 2D"
        self.radio_2d = QRadioButton(radio_2d_text)
        self.radio_3d.setChecked(True)
        self.radio_3d.toggled.connect(self.on_mode_changed)
        self.radio_2d.toggled.connect(self.on_mode_changed)
        
        input_layout.addWidget(self.line_edit)
        input_layout.addWidget(self.button)
        input_layout.addWidget(self.radio_3d)
        input_layout.addWidget(self.radio_2d)
        
        # повідомлення про помилку
        self.error_label = QLabel()
        self.error_label.setStyleSheet("color: red;")
        self.error_label.hide()
        
 
        self.visualization_frame = QFrame()
        self.visualization_layout = QVBoxLayout(self.visualization_frame)
        self.visualization_frame.hide()
        
        # таблиця
        self.create_molecule_table()
        
    
        self.create_action_buttons()
        
        # додаємо всі елементи
        self.main_tab_layout.addLayout(input_layout)
        self.main_tab_layout.addWidget(self.error_label, alignment=Qt.AlignTop)
        self.main_tab_layout.addWidget(self.info_label)
        self.main_tab_layout.addWidget(self.visualization_frame)
        self.main_tab_layout.addWidget(self.mol_table)
        self.main_tab_layout.addWidget(self.buttons_widget)

    def create_molecule_table(self):
        self.mol_table = QTableWidget()
        self.mol_table.setRowCount(20) 
        self.mol_table.setColumnCount(4)
        if self.lang == 'en':
            h_headers = ["Property", "Molecule 1", "Molecule 2", "Molecule 3"]
        elif self.lang == 'de':
            h_headers = ["Eigenschaft", "Molekül 1", "Molekül 2", "Molekül 3"]
        elif self.lang == 'uk':
            h_headers = ["Властивість", "Молекула 1", "Молекула 2", "Молекула 3"]
        elif self.lang == 'fr':
            h_headers = ["Propriété", "Molécule 1", "Molécule 2", "Molécule 3"]
        self.mol_table.setHorizontalHeaderLabels(h_headers)
        if self.lang == 'en':
            v_headers = [
                "Chemical formula", "Molecular weight", "LogP", "Number of atoms", "Number of bonds",
                "Aromatic atoms", "Rings", "H donors", "H acceptors",
                "E-state", "Sum of weights", "Surface area", "Number of C atoms",
                "Number of O atoms", "Number of N atoms", "Number of H atoms",
                "Number of halogens", "Polar surface area", "Allotropic logP",
                "Number of rotatable bonds"
            ]
        elif self.lang == 'de':
            v_headers = [
                "Chemische Formel", "Molekulargewicht", "LogP", "Anzahl der Atome", "Anzahl der Bindungen",
                "Aromatische Atome", "Ringe", "H-Spender", "H-Akzeptoren",
                "E-State", "Summe der Gewichte", "Oberflächenfläche", "Anzahl der C-Atome",
                "Anzahl der O-Atome", "Anzahl der N-Atome", "Anzahl der H-Atome",
                "Anzahl der Halogene", "Polare Oberflächenfläche", "Allotroper LogP",
                "Anzahl drehbarer Bindungen"
            ]
        elif self.lang == 'uk':
            v_headers = [
                "Хімічна формула", "Молекулярна вага", "LogP", "Число атомів", "Число зв'язків",
                "Ароматичні атоми", "Цикли", "Донори H", "Акцептори H",
                "E-state", "Сума ваг", "Площа поверхні", "Кількість атомів C",
                "Кількість атомів O", "Кількість атомів N", "Кількість атомів H",
                "Кількість галогенів", "Полярна площа поверхні", "Аллотропний logP",
                "Кількість обертових зв'язків"
            ]
        elif self.lang == 'fr':
            v_headers = [
                "Formule chimique", "Poids moléculaire", "LogP", "Nombre d'atomes", "Nombre de liaisons",
                "Atomes aromatiques", "Cycles", "Donneurs H", "Accepteurs H",
                "E-state", "Somme des poids", "Aire de surface", "Nombre d'atomes C",
                "Nombre d'atomes O", "Nombre d'atomes N", "Nombre d'atomes H",
                "Nombre d'halogènes", "Aire de surface polaire", "LogP allotrope",
                "Nombre de liaisons rotatives"
            ]
        self.mol_table.setVerticalHeaderLabels(v_headers)
        self.mol_table.hide()

    def create_action_buttons(self):
        self.line_edit_add = QLineEdit()
        if self.lang == 'en':
            add_placeholder = "Enter SMILES of second molecule..."
        elif self.lang == 'de':
            add_placeholder = "SMILES des zweiten Moleküls eingeben..."
        elif self.lang == 'uk':
            add_placeholder = "Введіть SMILES другої молекули..."
        elif self.lang == 'fr':
            add_placeholder = "Entrez SMILES de la deuxième molécule..."
        self.line_edit_add.setPlaceholderText(add_placeholder)
        if self.lang == 'en':
            add_button_text = "Add second molecule"
        elif self.lang == 'de':
            add_button_text = "Zweites Molekül hinzufügen"
        elif self.lang == 'uk':
            add_button_text = "Додати другу молекулу"
        elif self.lang == 'fr':
            add_button_text = "Ajouter deuxième molécule"
        self.button_add = QPushButton(add_button_text)
        self.button_add.clicked.connect(self.add_molecule)
        
        self.line_edit_add_3 = QLineEdit()
        if self.lang == 'en':
            add3_placeholder = "Enter SMILES of third molecule..."
        elif self.lang == 'de':
            add3_placeholder = "SMILES des dritten Moleküls eingeben..."
        elif self.lang == 'uk':
            add3_placeholder = "Введіть SMILES третьої молекули..."
        elif self.lang == 'fr':
            add3_placeholder = "Entrez SMILES de la troisième molécule..."
        self.line_edit_add_3.setPlaceholderText(add3_placeholder)
        if self.lang == 'en':
            add3_button_text = "Add third molecule"
        elif self.lang == 'de':
            add3_button_text = "Drittes Molekül hinzufügen"
        elif self.lang == 'uk':
            add3_button_text = "Додати третю молекулу"
        elif self.lang == 'fr':
            add3_button_text = "Ajouter troisième molécule"
        self.button_add_3 = QPushButton(add3_button_text)
        self.button_add_3.clicked.connect(self.add_molecule_3)
        
        if self.lang == 'en':
            export_text = "Export"
        elif self.lang == 'de':
            export_text = "Exportieren"
        elif self.lang == 'uk':
            export_text = "Експорт"
        elif self.lang == 'fr':
            export_text = "Exporter"
        self.export_button = QPushButton(export_text)
        self.export_button.clicked.connect(self.parent.export_molecule_advanced)
        
        if self.lang == 'en':
            fullscreen_3d_text = "Fullscreen 3D"
        elif self.lang == 'de':
            fullscreen_3d_text = "Vollbild 3D"
        elif self.lang == 'uk':
            fullscreen_3d_text = "Повноекранний 3D"
        elif self.lang == 'fr':
            fullscreen_3d_text = "Plein écran 3D"
        self.open_3d_button = QPushButton(fullscreen_3d_text)
        self.open_3d_button.clicked.connect(self.open_fullscreen_3d)
        
        if self.lang == 'en':
            fullscreen_table_text = "Fullscreen table"
        elif self.lang == 'de':
            fullscreen_table_text = "Vollbild Tabelle"
        elif self.lang == 'uk':
            fullscreen_table_text = "Повноекранна таблиця"
        elif self.lang == 'fr':
            fullscreen_table_text = "Tableau plein écran"
        self.open_table_button = QPushButton(fullscreen_table_text)
        self.open_table_button.clicked.connect(self.open_fullscreen_table)
        
        if self.lang == 'en':
            diagram_text = "Diagram"
        elif self.lang == 'de':
            diagram_text = "Diagramm"
        elif self.lang == 'uk':
            diagram_text = "Діаграма"
        elif self.lang == 'fr':
            diagram_text = "Diagramme"
        self.diagram_button = QPushButton(diagram_text)
        self.diagram_button.clicked.connect(self.show_comparison_diagram)
        
        if self.lang == 'en':
            clear_text = "Clear all"
        elif self.lang == 'de':
            clear_text = "Alles löschen"
        elif self.lang == 'uk':
            clear_text = "Очистити все"
        elif self.lang == 'fr':
            clear_text = "Effacer tout"
        self.clear_button = QPushButton(clear_text)
        self.clear_button.clicked.connect(self.clear_all_molecules)
        
        buttons_grid = QGridLayout()
        buttons_grid.addWidget(self.line_edit_add, 0, 0)
        buttons_grid.addWidget(self.button_add, 0, 1)
        buttons_grid.addWidget(self.line_edit_add_3, 0, 2)
        buttons_grid.addWidget(self.button_add_3, 0, 3)
        buttons_grid.addWidget(self.open_3d_button, 1, 0)
        buttons_grid.addWidget(self.export_button, 1, 1)
        buttons_grid.addWidget(self.open_table_button, 1, 2)
        buttons_grid.addWidget(self.diagram_button, 1, 3)
        buttons_grid.addWidget(self.clear_button, 1, 4)
        
        self.buttons_widget = QWidget()
        self.buttons_widget.setLayout(buttons_grid)
        self.buttons_widget.hide()

    def show_comparison_diagram(self):
        molecules_data = []
        
        for mol in [self.current_mol, self.current_mol_2, self.current_mol_3]:
            if mol:
                lang = self.lang
                if lang == 'en':
                    mol_data = {
                        'Molecular weight': Descriptors.MolWt(mol),
                        'LogP': Descriptors.MolLogP(mol),
                        'Number of atoms': mol.GetNumAtoms(),
                        'Number of bonds': mol.GetNumBonds(),
                        'Aromatic atoms': len([a for a in mol.GetAtoms() if a.GetIsAromatic()]),
                        'Rings': Chem.GetSSSR(mol),
                        'H donors': Descriptors.NumHDonors(mol),
                        'H acceptors': Descriptors.NumHAcceptors(mol),
                    }
                elif lang == 'de':
                    mol_data = {
                        'Molekulargewicht': Descriptors.MolWt(mol),
                        'LogP': Descriptors.MolLogP(mol),
                        'Anzahl der Atome': mol.GetNumAtoms(),
                        'Anzahl der Bindungen': mol.GetNumBonds(),
                        'Aromatische Atome': len([a for a in mol.GetAtoms() if a.GetIsAromatic()]),
                        'Ringe': Chem.GetSSSR(mol),
                        'H-Spender': Descriptors.NumHDonors(mol),
                        'H-Akzeptoren': Descriptors.NumHAcceptors(mol),
                    }
                elif lang == 'uk':
                    mol_data = {
                        'Молекулярна вага': Descriptors.MolWt(mol),
                        'LogP': Descriptors.MolLogP(mol),
                        'Число атомів': mol.GetNumAtoms(),
                        'Число зв\'язків': mol.GetNumBonds(),
                        'Ароматичні атоми': len([a for a in mol.GetAtoms() if a.GetIsAromatic()]),
                        'Цикли': Chem.GetSSSR(mol),
                        'Донори H': Descriptors.NumHDonors(mol),
                        'Акцептори H': Descriptors.NumHAcceptors(mol),
                    }
                elif lang == 'fr':
                    mol_data = {
                        'Poids moléculaire': Descriptors.MolWt(mol),
                        'LogP': Descriptors.MolLogP(mol),
                        'Nombre d\'atomes': mol.GetNumAtoms(),
                        'Nombre de liaisons': mol.GetNumBonds(),
                        'Atomes aromatiques': len([a for a in mol.GetAtoms() if a.GetIsAromatic()]),
                        'Cycles': Chem.GetSSSR(mol),
                        'Donneurs H': Descriptors.NumHDonors(mol),
                        'Accepteurs H': Descriptors.NumHAcceptors(mol),
                    }
                molecules_data.append(mol_data)
        
        if len(molecules_data) < 2:
            lang = self.lang
            if lang == 'en':
                warning_text = "To create a diagram, at least two molecules are needed"
            elif lang == 'de':
                warning_text = "Um ein Diagramm zu erstellen, sind mindestens zwei Moleküle erforderlich"
            elif lang == 'uk':
                warning_text = "Для створення діаграми потрібно щонайменше дві молекули"
            elif lang == 'fr':
                warning_text = "Pour créer un diagramme, au moins deux molécules sont nécessaires"
            warning_title = "Warning" if lang == 'en' else "Warnung" if lang == 'de' else "Попередження" if lang == 'uk' else "Avertissement" if lang == 'fr' else "Warning"
            QMessageBox.warning(self, warning_title, warning_text)
            return
            
        diagram_tab = ComparisonDiagramTab(molecules_data, self.parent)
        lang = self.lang
        if lang == 'en':
            tab_title = "Comparative Diagram"
        elif lang == 'de':
            tab_title = "Vergleichsdiagramm"
        elif lang == 'uk':
            tab_title = "Порівняльна діаграма"
        elif lang == 'fr':
            tab_title = "Diagramme comparatif"
        self.parent.tabs.addTab(diagram_tab, tab_title)
        self.parent.tabs.setCurrentIndex(self.parent.tabs.count()-1)

    def clear_visualization(self):
        while self.visualization_layout.count():
            item = self.visualization_layout.takeAt(0)
            widget = item.widget()
            if widget:
                widget.deleteLater()
            else:
                nested_layout = item.layout()
                if nested_layout:
                    while nested_layout.count():
                        nested_item = nested_layout.takeAt(0)
                        nested_widget = nested_item.widget()
                        if nested_widget:
                            nested_widget.deleteLater()
                    nested_layout.deleteLater()

    def is_valid_smiles(self, smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False
            Chem.SanitizeMol(mol, Chem.SANITIZE_PROPERTIES)
            return True
        except:
            return False

    def check_smiles(self):
        smiles = self.line_edit.text().strip()
        
        if not smiles:
            lang = self.lang
            if lang == 'en':
                error_msg = "Enter SMILES!"
            elif lang == 'de':
                error_msg = "SMILES eingeben!"
            elif lang == 'uk':
                error_msg = "Введіть SMILES!"
            elif lang == 'fr':
                error_msg = "Entrez SMILES!"
            self.show_error(error_msg)
            return
            
        if not self.is_valid_smiles(smiles):
            lang = self.lang
            if lang == 'en':
                error_msg = "Invalid SMILES!"
            elif lang == 'de':
                error_msg = "Ungültiges SMILES!"
            elif lang == 'uk':
                error_msg = "Невірний SMILES!"
            elif lang == 'fr':
                error_msg = "SMILES invalide!"
            self.show_error(error_msg)
            return
            
        self.error_label.hide()
        
        try:
            self.current_mol = Chem.MolFromSmiles(smiles)
            if not self.current_mol:
                lang = self.lang
                if lang == 'en':
                    error_msg = "Error creating molecule"
                elif lang == 'de':
                    error_msg = "Fehler beim Erstellen des Moleküls"
                elif lang == 'uk':
                    error_msg = "Помилка створення молекули"
                elif lang == 'fr':
                    error_msg = "Erreur lors de la création de la molécule"
                self.show_error(error_msg)
                return
                
            self.current_mol = Chem.AddHs(self.current_mol)
            AllChem.EmbedMolecule(self.current_mol)
            AllChem.UFFOptimizeMolecule(self.current_mol)
            self.current_mol_block = Chem.MolToMolBlock(self.current_mol)
            
            self.update_molecule_info()
            
            if self.radio_3d.isChecked():
                self.show_3d_view()
            else:
                self.show_2d_view()
                
            # спрятати інфо текст
            self.info_label.hide()
            self.visualization_frame.show()
        except Exception as e:
            lang = self.lang
            if lang == 'en':
                error_msg = "Error processing molecule: "
            elif lang == 'de':
                error_msg = "Fehler bei der Verarbeitung des Moleküls: "
            elif lang == 'uk':
                error_msg = "Помилка обробки молекули: "
            elif lang == 'fr':
                error_msg = "Erreur lors du traitement de la molécule : "
            self.show_error(f"{error_msg}{str(e)}")

    def update_molecule_info(self):
        molecules = [self.current_mol, self.current_mol_2, self.current_mol_3]
        props = []
        
        for mol in molecules:
            if mol:
                formula = rdMolDescriptors.CalcMolFormula(mol)
                
                # Additional properties
                num_carbon = len([a for a in mol.GetAtoms() if a.GetSymbol() == 'C'])
                num_oxygen = len([a for a in mol.GetAtoms() if a.GetSymbol() == 'O'])
                num_nitrogen = len([a for a in mol.GetAtoms() if a.GetSymbol() == 'N'])
                num_hydrogen = len([a for a in mol.GetAtoms() if a.GetSymbol() == 'H'])
                num_halogens = len([a for a in mol.GetAtoms() if a.GetSymbol() in ['F', 'Cl', 'Br', 'I']])
                
                mol_props = [
                    formula,
                    f"{Descriptors.MolWt(mol):.2f} g/mol",
                    f"{Descriptors.MolLogP(mol):.2f}",
                    str(mol.GetNumAtoms()),
                    str(mol.GetNumBonds()),
                    str(len([a for a in mol.GetAtoms() if a.GetIsAromatic()])),
                    str(Chem.GetSSSR(mol)),
                    str(Descriptors.NumHDonors(mol)),
                    str(Descriptors.NumHAcceptors(mol)),
                    f"{Descriptors.EState_VSA1(mol):.2f}",
                    f"{sum(a.GetMass() for a in mol.GetAtoms()):.2f} g",
                    f"{rdMolDescriptors.CalcTPSA(mol):.2f} Å²",
                    str(num_carbon),
                    str(num_oxygen),
                    str(num_nitrogen),
                    str(num_hydrogen),
                    str(num_halogens),
                    f"{rdMolDescriptors.CalcTPSA(mol):.2f} Å²",  
                    f"{Descriptors.MolLogP(mol):.2f}",  
                    str(Descriptors.NumRotatableBonds(mol))  
                ]
            else:
                mol_props = [""] * 20
            
            props.append(mol_props)
        
        for row in range(20):
            for col in range(3):
                self.mol_table.setItem(row, col+1, QTableWidgetItem(props[col][row]))
        
        self.mol_table.resizeColumnsToContents()
        if any(molecules):
            self.mol_table.show()
    # основні функції програми
    def add_molecule(self):
        smiles = self.line_edit_add.text().strip()
        
        if not smiles:
            lang = self.lang
            if lang == 'en':
                error_msg = "Enter SMILES for second molecule!"
            elif lang == 'de':
                error_msg = "SMILES für zweites Molekül eingeben!"
            elif lang == 'uk':
                error_msg = "Введіть SMILES для другої молекули!"
            elif lang == 'fr':
                error_msg = "Entrez SMILES pour la deuxième molécule!"
            self.show_error(error_msg)
            return
            
        if not self.is_valid_smiles(smiles):
            lang = self.lang
            if lang == 'en':
                error_msg = "Invalid SMILES for second molecule!"
            elif lang == 'de':
                error_msg = "Ungültiges SMILES für zweites Molekül!"
            elif lang == 'uk':
                error_msg = "Невірний SMILES для другої молекули!"
            elif lang == 'fr':
                error_msg = "SMILES invalide pour la deuxième molécule!"
            self.show_error(error_msg)
            return
            
        self.error_label.hide()
        
        try:
            self.current_mol_2 = Chem.MolFromSmiles(smiles)
            if not self.current_mol_2:
                lang = self.lang
                if lang == 'en':
                    error_msg = "Error creating second molecule"
                elif lang == 'de':
                    error_msg = "Fehler beim Erstellen des zweiten Moleküls"
                elif lang == 'uk':
                    error_msg = "Помилка створення другої молекули"
                elif lang == 'fr':
                    error_msg = "Erreur lors de la création de la deuxième molécule"
                self.show_error(error_msg)
                return
                
            self.current_mol_2 = Chem.AddHs(self.current_mol_2)
            AllChem.EmbedMolecule(self.current_mol_2)
            AllChem.UFFOptimizeMolecule(self.current_mol_2)
            
            self.update_molecule_info()
            
            if self.radio_3d.isChecked():
                self.show_3d_view()
            else:
                self.show_2d_view()
                
            self.buttons_widget.show()
        except Exception as e:
            lang = self.lang
            if lang == 'en':
                error_msg = "Error processing second molecule: "
            elif lang == 'de':
                error_msg = "Fehler bei der Verarbeitung des zweiten Moleküls: "
            elif lang == 'uk':
                error_msg = "Помилка обробки другої молекули: "
            elif lang == 'fr':
                error_msg = "Erreur lors du traitement de la deuxième molécule : "
            self.show_error(f"{error_msg}{str(e)}")

    def add_molecule_3(self):
        smiles = self.line_edit_add_3.text().strip()
        
        if not smiles:
            lang = self.lang
            if lang == 'en':
                error_msg = "Enter SMILES for third molecule!"
            elif lang == 'de':
                error_msg = "SMILES für drittes Molekül eingeben!"
            elif lang == 'uk':
                error_msg = "Введіть SMILES для третьої молекули!"
            elif lang == 'fr':
                error_msg = "Entrez SMILES pour la troisième molécule!"
            self.show_error(error_msg)
            return
            
        if not self.is_valid_smiles(smiles):
            lang = self.lang
            if lang == 'en':
                error_msg = "Invalid SMILES for third molecule!"
            elif lang == 'de':
                error_msg = "Ungültiges SMILES für drittes Molekül!"
            elif lang == 'uk':
                error_msg = "Невірний SMILES для третьої молекули!"
            elif lang == 'fr':
                error_msg = "SMILES invalide pour la troisième molécule!"
            self.show_error(error_msg)
            return
            
        self.error_label.hide()
        
        try:
            self.current_mol_3 = Chem.MolFromSmiles(smiles)
            if not self.current_mol_3:
                lang = self.lang
                if lang == 'en':
                    error_msg = "Error creating third molecule"
                elif lang == 'de':
                    error_msg = "Fehler beim Erstellen des dritten Moleküls"
                elif lang == 'uk':
                    error_msg = "Помилка створення третьої молекули"
                elif lang == 'fr':
                    error_msg = "Erreur lors de la création de la troisième molécule"
                self.show_error(error_msg)
                return
                
            self.current_mol_3 = Chem.AddHs(self.current_mol_3)
            AllChem.EmbedMolecule(self.current_mol_3)
            AllChem.UFFOptimizeMolecule(self.current_mol_3)
            
            self.update_molecule_info()
            
            if self.radio_3d.isChecked():
                self.show_3d_view()
            else:
                self.show_2d_view()
                
            self.buttons_widget.show()
        except Exception as e:
            lang = self.lang
            if lang == 'en':
                error_msg = "Error processing third molecule: "
            elif lang == 'de':
                error_msg = "Fehler bei der Verarbeitung des dritten Moleküls: "
            elif lang == 'uk':
                error_msg = "Помилка обробки третьої молекули: "
            elif lang == 'fr':
                error_msg = "Erreur lors du traitement de la troisième molécule : "
            self.show_error(f"{error_msg}{str(e)}")

    def clear_all_molecules(self):
        self.current_mol = None
        self.current_mol_2 = None
        self.current_mol_3 = None
        
        self.line_edit.clear()
        self.line_edit_add.clear()
        self.line_edit_add_3.clear()
        
        self.clear_visualization()
        self.update_molecule_info()
        self.mol_table.hide()
        self.buttons_widget.hide()
        self.visualization_frame.hide()
        self.info_label.show()

    def show_3d_view(self):
        self.clear_visualization()
        horizontal_layout = QHBoxLayout()

        molecules = [self.current_mol, self.current_mol_2, self.current_mol_3]
        
        for i, mol in enumerate(molecules):
            if mol:
                try:
                    view = QWebEngineView()
                    view.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
                    view.setMinimumSize(300, 300)

                    viewer = py3Dmol.view(width=300, height=300)
                    mol_block = Chem.MolToMolBlock(mol)
                    viewer.addModel(mol_block, f"mol_{i}")
                    
                    if self.show_hydrogens:
                        viewer.setStyle({self.visualization_style: {}})
                    else:
                        viewer.setStyle({self.visualization_style: {'sel': {'elem': {'H': False}}}})
                    
                    viewer.zoomTo()
                
                    view.setHtml(viewer._make_html())
                    horizontal_layout.addWidget(view)
                except Exception as e:
                    lang = self.lang
                    if lang == 'en':
                        error_msg = f"Error visualizing molecule {i+1}: "
                    elif lang == 'de':
                        error_msg = f"Fehler bei der Visualisierung des Moleküls {i+1}: "
                    elif lang == 'uk':
                        error_msg = f"Помилка візуалізації молекули {i+1}: "
                    elif lang == 'fr':
                        error_msg = f"Erreur de visualisation de la molécule {i+1} : "
                    self.show_error(f"{error_msg}{str(e)}")

        self.visualization_layout.addLayout(horizontal_layout)
        self.buttons_widget.show()

    def show_2d_view(self):
        self.clear_visualization()
        horizontal_layout = QHBoxLayout()
    
        molecules = [self.current_mol, self.current_mol_2, self.current_mol_3]
        
        for i, mol in enumerate(molecules):
            if mol:
                try:
                    pixmap = self.mol_to_qpixmap(Chem.MolToSmiles(mol))
                    if pixmap:
                        label_2d = QLabel()
                        label_2d.setFixedSize(300, 300)
                        label_2d.setPixmap(pixmap)
                        horizontal_layout.addWidget(label_2d)
                except Exception as e:
                    lang = self.lang
                    if lang == 'en':
                        error_msg = f"Error creating 2D image of molecule {i+1}: "
                    elif lang == 'de':
                        error_msg = f"Fehler beim Erstellen des 2D-Bildes des Moleküls {i+1}: "
                    elif lang == 'uk':
                        error_msg = f"Помилка створення 2D зображення молекули {i+1}: "
                    elif lang == 'fr':
                        error_msg = f"Erreur lors de la création de l'image 2D de la molécule {i+1} : "
                    self.show_error(f"{error_msg}{str(e)}")
    
        self.visualization_layout.addLayout(horizontal_layout)
        self.buttons_widget.show()

    def mol_to_qpixmap(self, smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return None
            
            img = Draw.MolToImage(mol, size=(300, 300))
            buffer = BytesIO()
            img.save(buffer, format="PNG")
            buffer.seek(0)
        
            pixmap = QPixmap()
            pixmap.loadFromData(buffer.getvalue())
            return pixmap
        except:
            return None

    def open_fullscreen_3d(self):
        new_tab = QWidget()
        layout = QVBoxLayout(new_tab)
        
        # повний екран для 3д
        button_layout = QHBoxLayout()
        
        # зберегти усі три моделі
        lang = self.lang
        if lang == 'en':
            save_all_text = "Save all 3 models"
        elif lang == 'de':
            save_all_text = "Alle 3 Modelle speichern"
        elif lang == 'uk':
            save_all_text = "Зберегти всі 3 моделі"
        elif lang == 'fr':
            save_all_text = "Sauvegarder tous les 3 modèles"
        save_all_button = QPushButton(save_all_text)
        save_all_button.clicked.connect(lambda: self.save_all_3_models())
        button_layout.addWidget(save_all_button)
        
        layout.addLayout(button_layout)
        
        # дадати 3д
        view_layout = QHBoxLayout()
        
        molecules = [self.current_mol, self.current_mol_2, self.current_mol_3]
        
        for i, mol in enumerate(molecules):
            if mol:
                try:
                    view = QWebEngineView()
                    view.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
                    
                    viewer = py3Dmol.view(width=600, height=800)
                    mol_block = Chem.MolToMolBlock(mol)
                    viewer.addModel(mol_block, f"mol_{i}")
                    
                    if self.show_hydrogens:
                        viewer.setStyle({self.visualization_style: {}})
                    else:
                        viewer.setStyle({self.visualization_style: {'sel': {'elem': {'H': False}}}})
                        
                    viewer.zoomTo()
                    
                    view.setHtml(viewer._make_html())
                    view_layout.addWidget(view)
                except Exception as e:
                    lang = self.lang
                    if lang == 'en':
                        error_msg = f"Error in fullscreen visualization of molecule {i+1}: "
                    elif lang == 'de':
                        error_msg = f"Fehler bei der Vollbildvisualisierung des Moleküls {i+1}: "
                    elif lang == 'uk':
                        error_msg = f"Помилка повноекранної візуалізації молекули {i+1}: "
                    elif lang == 'fr':
                        error_msg = f"Erreur dans la visualisation plein écran de la molécule {i+1} : "
                    self.show_error(f"{error_msg}{str(e)}")
        
        layout.addLayout(view_layout)
        
        lang = self.lang
        if lang == 'en':
            tab_title = "3D Fullscreen"
        elif lang == 'de':
            tab_title = "3D Vollbild"
        elif lang == 'uk':
            tab_title = "3D Повний екран"
        elif lang == 'fr':
            tab_title = "3D Plein écran"
        self.parent.tabs.addTab(new_tab, tab_title)
        self.parent.tabs.setCurrentIndex(self.parent.tabs.count()-1)

    def open_fullscreen_table(self):
        new_tab = QWidget()
        layout = QVBoxLayout(new_tab)
        
        # полний екран для таблиць
        button_layout = QHBoxLayout()
        
        # кнопка перевірити схожість
        lang = self.lang
        if lang == 'en':
            check_columns_text = "Check columns"
        elif lang == 'de':
            check_columns_text = "Spalten überprüfen"
        elif lang == 'uk':
            check_columns_text = "Перевірити стовпчики"
        elif lang == 'fr':
            check_columns_text = "Vérifier les colonnes"
        check_columns_button = QPushButton(check_columns_text)
        check_columns_button.clicked.connect(lambda: self.check_table_columns(fullscreen_table))
        button_layout.addWidget(check_columns_button)
        
        # експортувати таблицю
        if lang == 'en':
            export_table_text = "Export table"
        elif lang == 'de':
            export_table_text = "Tabelle exportieren"
        elif lang == 'uk':
            export_table_text = "Експорт таблиці"
        elif lang == 'fr':
            export_table_text = "Exporter le tableau"
        export_table_button = QPushButton(export_table_text)
        export_table_button.clicked.connect(lambda: self.export_table_to_file(fullscreen_table))
        button_layout.addWidget(export_table_button)
        
        layout.addLayout(button_layout)
        
        fullscreen_table = QTableWidget()
        fullscreen_table.setRowCount(20)
        fullscreen_table.setColumnCount(4)
        if lang == 'en':
            h_headers = ["Property", "Molecule 1", "Molecule 2", "Molecule 3"]
        elif lang == 'de':
            h_headers = ["Eigenschaft", "Molekül 1", "Molekül 2", "Molekül 3"]
        elif lang == 'uk':
            h_headers = ["Властивість", "Молекула 1", "Молекула 2", "Молекула 3"]
        elif lang == 'fr':
            h_headers = ["Propriété", "Molécule 1", "Molécule 2", "Molécule 3"]
        fullscreen_table.setHorizontalHeaderLabels(h_headers)
        if lang == 'en':
            v_headers = [
                "Chemical formula", "Molecular weight", "LogP", "Number of atoms", "Number of bonds",
                "Aromatic atoms", "Rings", "H donors", "H acceptors",
                "E-state", "Sum of weights", "Surface area", "Number of C atoms",
                "Number of O atoms", "Number of N atoms", "Number of H atoms",
                "Number of halogens", "Polar surface area", "Allotropic logP",
                "Number of rotatable bonds"
            ]
        elif lang == 'de':
            v_headers = [
                "Chemische Formel", "Molekulargewicht", "LogP", "Anzahl der Atome", "Anzahl der Bindungen",
                "Aromatische Atome", "Ringe", "H-Spender", "H-Akzeptoren",
                "E-State", "Summe der Gewichte", "Oberflächenfläche", "Anzahl der C-Atome",
                "Anzahl der O-Atome", "Anzahl der N-Atome", "Anzahl der H-Atome",
                "Anzahl der Halogene", "Polare Oberflächenfläche", "Allotroper LogP",
                "Anzahl drehbarer Bindungen"
            ]
        elif lang == 'uk':
            v_headers = [
                "Хімічна формула", "Молекулярна вага", "LogP", "Число атомів", "Число зв'язків",
                "Ароматичні атоми", "Цикли", "Донори H", "Акцептори H",
                "E-state", "Сума ваг", "Площа поверхні", "Кількість атомів C",
                "Кількість атомів O", "Кількість атомів N", "Кількість атомів H",
                "Кількість галогенів", "Полярна площа поверхні", "Аллотропний logP",
                "Кількість обертових зв'язків"
            ]
        elif lang == 'fr':
            v_headers = [
                "Formule chimique", "Poids moléculaire", "LogP", "Nombre d'atomes", "Nombre de liaisons",
                "Atomes aromatiques", "Cycles", "Donneurs H", "Accepteurs H",
                "E-state", "Somme des poids", "Aire de surface", "Nombre d'atomes C",
                "Nombre d'atomes O", "Nombre d'atomes N", "Nombre d'atomes H",
                "Nombre d'halogènes", "Aire de surface polaire", "LogP allotrope",
                "Nombre de liaisons rotatives"
            ]
        fullscreen_table.setVerticalHeaderLabels(v_headers)
        

        for row in range(20):
            for col in range(4):
                item = self.mol_table.item(row, col)
                if item:
                    fullscreen_table.setItem(row, col, QTableWidgetItem(item.text()))
        
        fullscreen_table.resizeColumnsToContents()
        layout.addWidget(fullscreen_table)
        
        lang = self.lang
        if lang == 'en':
            tab_title = "Fullscreen Table"
        elif lang == 'de':
            tab_title = "Vollbild Tabelle"
        elif lang == 'uk':
            tab_title = "Повноекранна таблиця"
        elif lang == 'fr':
            tab_title = "Tableau plein écran"
        self.parent.tabs.addTab(new_tab, tab_title)
        self.parent.tabs.setCurrentIndex(self.parent.tabs.count()-1)

    def on_mode_changed(self):
        if self.current_mol or self.current_mol_2 or self.current_mol_3:
            if self.radio_3d.isChecked():
                self.show_3d_view()
            else:
                self.show_2d_view()

    def show_error(self, message):
        self.error_label.setText(message)
        self.error_label.show()

    def set_visualization_style(self, style):
        self.visualization_style = style
        if self.radio_3d.isChecked() and (self.current_mol or self.current_mol_2 or self.current_mol_3):
            self.show_3d_view()

    def toggle_hydrogens(self):
        self.show_hydrogens = not self.show_hydrogens
        if self.radio_3d.isChecked() and (self.current_mol or self.current_mol_2 or self.current_mol_3):
            self.show_3d_view()

    def get_data(self):
        data = {
            "smiles": self.line_edit.text(),
            "smiles2": self.line_edit_add.text() if self.line_edit_add else "",
            "smiles3": self.line_edit_add_3.text() if self.line_edit_add_3 else "",
        }
        return data

    def load_data(self, data):
        self.line_edit.setText(data.get("smiles", ""))
        if "smiles2" in data:
            self.line_edit_add.setText(data["smiles2"])
            self.add_molecule()
        if "smiles3" in data:
            self.line_edit_add_3.setText(data["smiles3"])
            self.add_molecule_3()


    def save_all_3_models(self):
        try:
            molecules = [self.current_mol, self.current_mol_2, self.current_mol_3]
            if not any(molecules):
                lang = self.lang
                warning_title = "Warning" if lang == 'en' else "Warnung" if lang == 'de' else "Попередження" if lang == 'uk' else "Avertissement" if lang == 'fr' else "Warning"
                warning_msg = "No molecules to save" if lang == 'en' else "Keine Moleküle zum Speichern" if lang == 'de' else "Немає молекул для збереження" if lang == 'uk' else "Aucune molécule à sauvegarder" if lang == 'fr' else "No molecules to save"
                QMessageBox.warning(self, warning_title, warning_msg)
                return
                
            lang = self.lang
            title = "Save all 3 models" if lang == 'en' else "Alle 3 Modelle speichern" if lang == 'de' else "Зберегти всі 3 моделі" if lang == 'uk' else "Sauvegarder tous les 3 modèles" if lang == 'fr' else "Save all 3 models"
            filename, _ = QFileDialog.getSaveFileName(self, title, "", "SDF Files (*.sdf);;MOL Files (*.mol);;All Files (*)")
            if filename:
                if filename.endswith('.sdf') or '.sdf' in filename:
                    writer = Chem.SDWriter(filename)
                    for i, mol in enumerate(molecules):
                        if mol:
                            writer.write(mol)
                    writer.close()
                else:
                    with open(filename, 'w') as f:
                        for i, mol in enumerate(molecules):
                            if mol:
                                mol_block = Chem.MolToMolBlock(mol)
                                f.write(mol_block)
                                f.write("\n$$$$\n")
                lang = self.lang
                success_title = "Success" if lang == 'en' else "Erfolg" if lang == 'de' else "Успіх" if lang == 'uk' else "Succès" if lang == 'fr' else "Success"
                success_msg = "All 3 models successfully saved" if lang == 'en' else "Alle 3 Modelle erfolgreich gespeichert" if lang == 'de' else "Всі 3 моделі успішно збережено" if lang == 'uk' else "Tous les 3 modèles sauvegardés avec succès" if lang == 'fr' else "All 3 models successfully saved"
                QMessageBox.information(self, success_title, success_msg)
        except Exception as e:
            lang = self.lang
            error_title = "Error" if lang == 'en' else "Fehler" if lang == 'de' else "Помилка" if lang == 'uk' else "Erreur" if lang == 'fr' else "Error"
            error_msg = "Error saving: " if lang == 'en' else "Fehler beim Speichern: " if lang == 'de' else "Помилка збереження: " if lang == 'uk' else "Erreur lors de la sauvegarde : " if lang == 'fr' else "Error saving: "
            QMessageBox.warning(self, error_title, f"{error_msg}{str(e)}")

    def check_table_columns(self, table):
        try:
            for row in range(table.rowCount()):
                values = []
                for col in range(1, table.columnCount()): 
                    item = table.item(row, col)
                    if item and item.text():
                        values.append(item.text())
                
                if len(values) > 1 and all(v == values[0] for v in values):
                    for col in range(1, table.columnCount()):
                        item = table.item(row, col)
                        if item:
                            item.setBackground(Qt.green)
                else:
                    for col in range(1, table.columnCount()):
                        item = table.item(row, col)
                        if item:
                            item.setBackground(Qt.white)
        except Exception as e:
            lang = self.lang
            error_title = "Error" if lang == 'en' else "Fehler" if lang == 'de' else "Помилка" if lang == 'uk' else "Erreur" if lang == 'fr' else "Error"
            error_msg = "Error checking columns: " if lang == 'en' else "Fehler bei der Überprüfung der Spalten: " if lang == 'de' else "Помилка перевірки стовпчиків: " if lang == 'uk' else "Erreur lors de la vérification des colonnes : " if lang == 'fr' else "Error checking columns: "
            QMessageBox.warning(self, error_title, f"{error_msg}{str(e)}")

    def export_table_to_file(self, table):
        try:
            lang = self.lang
            title = "Export table" if lang == 'en' else "Tabelle exportieren" if lang == 'de' else "Експорт таблиці" if lang == 'uk' else "Exporter le tableau" if lang == 'fr' else "Export table"
            filename, _ = QFileDialog.getSaveFileName(self, title, "", "CSV Files (*.csv);;Excel Files (*.xlsx);;All Files (*)")
            if filename:
                if filename.endswith('.csv'):
                    self.export_table_to_csv(table, filename)
                elif filename.endswith('.xlsx'):
                    self.export_table_to_excel(table, filename)
                success_title = "Success" if lang == 'en' else "Erfolg" if lang == 'de' else "Успіх" if lang == 'uk' else "Succès" if lang == 'fr' else "Success"
                success_msg = "Table successfully exported" if lang == 'en' else "Tabelle erfolgreich exportiert" if lang == 'de' else "Таблицю успішно експортовано" if lang == 'uk' else "Tableau exporté avec succès" if lang == 'fr' else "Table successfully exported"
                QMessageBox.information(self, success_title, success_msg)
        except Exception as e:
            lang = self.lang
            error_title = "Error" if lang == 'en' else "Fehler" if lang == 'de' else "Помилка" if lang == 'uk' else "Erreur" if lang == 'fr' else "Error"
            error_msg = "Error exporting table: " if lang == 'en' else "Fehler beim Exportieren der Tabelle: " if lang == 'de' else "Помилка експорту таблиці: " if lang == 'uk' else "Erreur lors de l'export du tableau : " if lang == 'fr' else "Error exporting table: "
            QMessageBox.warning(self, error_title, f"{error_msg}{str(e)}")

    def export_table_to_csv(self, table, filename):
        try:
            with open(filename, 'w', encoding='utf-8') as f:

                headers = []
                for col in range(table.columnCount()):
                    headers.append(table.horizontalHeaderItem(col).text())
                f.write(','.join(headers) + '\n')
                
                for row in range(table.rowCount()):
                    row_data = []
                    for col in range(table.columnCount()):
                        item = table.item(row, col)
                        if item:
                            row_data.append(item.text())
                        else:
                            row_data.append('')
                    f.write(','.join(row_data) + '\n')
        except Exception as e:
            raise e

    def export_table_to_excel(self, table, filename):
        try:
            import pandas as pd
            data = []
            headers = []
            
            for col in range(table.columnCount()):
                headers.append(table.horizontalHeaderItem(col).text())
            
            for row in range(table.rowCount()):
                row_data = []
                for col in range(table.columnCount()):
                    item = table.item(row, col)
                    if item:
                        row_data.append(item.text())
                    else:
                        row_data.append('')
                data.append(row_data)
            df = pd.DataFrame(data, columns=headers)
            df.to_excel(filename, index=False)
        except Exception as e:
            raise e
