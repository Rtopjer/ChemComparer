# біблиотеки
from PyQt5.QtWidgets import (QApplication, QMainWindow, QVBoxLayout, QWidget, 
                            QLineEdit, QLabel, QPushButton, QFileDialog, 
                            QHBoxLayout, QRadioButton, QTableWidget, 
                            QTableWidgetItem, QTabWidget, QSizePolicy, QFrame, QMenuBar, QAction, QMessageBox, QMenu, QDialog, QScrollArea, QGridLayout, QStatusBar, QTabBar)
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtCore import Qt, QTranslator, QLocale, QSettings
from PyQt5.QtGui import QPixmap, QIcon
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors, rdMolDescriptors
from rdkit import RDLogger
import py3Dmol
from io import BytesIO
import sys
import json
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import pandas as pd
from dialogs import AboutDialog, SmilesHelpDialog
from widgets import ComparisonDiagramTab, MainTabWidget

RDLogger.DisableLog('rdApp.*')
#основне окно
class MoleculeViewer(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("ChemComparer")
        self.resize(1400, 900)
        self.setMinimumSize(1000, 700)
        
        self.setup_macos_icon()
        
        # налаштування 
        self.settings = QSettings("ChemComparer", "AppSettings")
        self.current_language = self.settings.value("language", "en", type=str)
        
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        
        self.create_menubar()
        
        self.main_layout = QVBoxLayout(self.central_widget)
        
        self.tabs = QTabWidget()
        self.tabs.setTabPosition(QTabWidget.West)
        self.tabs.setTabsClosable(True)
        self.tabs.tabCloseRequested.connect(self.close_tab)
        self.main_layout.addWidget(self.tabs)
        
        # начальна вкладка
        self.add_new_main_tab()
        self.tabs.tabBar().setTabButton(0, QTabBar.LeftSide, None)

        # прикольна штука
        self.setAcceptDrops(True)

        # горячі клафіши
        self.toggle_fullscreen_action = QAction(self)
        self.toggle_fullscreen_action.setShortcut("F11")
        self.toggle_fullscreen_action.triggered.connect(self.toggle_fullscreen)
        self.addAction(self.toggle_fullscreen_action)

        # 
        self.statusBar = QStatusBar()
        self.setStatusBar(self.statusBar)
        if self.current_language == 'en':
            status_msg = "ChemComparer 0.8.1 - To start, enter SMILES in the field above"
        elif self.current_language == 'de':
            status_msg = "ChemComparer 0.8.1 - Um zu starten, geben Sie SMILES in das Feld oben ein"
        elif self.current_language == 'uk':
            status_msg = "ChemComparer 0.8.1 - Для старту введіть SMILES у поле вище"
        elif self.current_language == 'fr':
            status_msg = "ChemComparer 0.8.1 - Pour commencer, entrez SMILES dans le champ ci-dessus"
        self.statusBar.showMessage(status_msg)
        
        # загруска теми
        self.load_theme_settings()

    def setup_macos_icon(self):
        icon_path = "icon.png"
        if os.path.exists(icon_path):
            self.setWindowIcon(QIcon(icon_path))
            if sys.platform == "darwin":
                app.setWindowIcon(QIcon(icon_path))

    def load_theme_settings(self):
        theme = self.settings.value("theme", "light", type=str)
        
        if theme == "dark":
            self.set_dark_theme()
        elif theme == "gray":
            self.set_gray_theme()
        else:
            self.set_light_theme()  # дефолтна тема

    def save_theme_settings(self, theme):
        """Saves theme settings"""
        self.settings.setValue("theme", theme)

    def add_new_main_tab(self):
        main_tab = MainTabWidget(self)
        lang = self.current_language
        if lang == 'en':
            tab_title = "Main"
        elif lang == 'de':
            tab_title = "Haupt"
        elif lang == 'uk':
            tab_title = "Головна"
        elif lang == 'fr':
            tab_title = "Principal"
        self.tabs.addTab(main_tab, tab_title)

    def close_tab(self, index):
        if self.tabs.count() > 1 and index > 0:  # незакриваємая вкладка
            self.tabs.removeTab(index)

    def close_all_tabs(self):
        while self.tabs.count() > 1:
            self.tabs.removeTab(1)

    def toggle_fullscreen(self):
        if self.isFullScreen():
            self.showNormal()
        else:
            self.showFullScreen()

    def create_menubar(self):
        menubar = self.menuBar()
        
        # файл меню
        lang = self.current_language
        file_menu_text = "File" if lang == 'en' else "Datei" if lang == 'de' else "Файл" if lang == 'uk' else "Fichier" if lang == 'fr' else "File"
        file_menu = menubar.addMenu(file_menu_text)
        
        new_tab_text = "New Tab" if lang == 'en' else "Neuer Tab" if lang == 'de' else "Нова вкладка" if lang == 'uk' else "Nouvel onglet" if lang == 'fr' else "New Tab"
        new_tab_action = QAction(new_tab_text, self)
        new_tab_action.triggered.connect(self.add_new_main_tab)
        file_menu.addAction(new_tab_action)

        save_text = "Save As" if lang == 'en' else "Speichern als" if lang == 'de' else "Зберегти як" if lang == 'uk' else "Enregistrer sous" if lang == 'fr' else "Save As"
        save_action = QAction(save_text, self)
        save_action.triggered.connect(self.save_data)
        save_action.setShortcut("Ctrl+S")
        file_menu.addAction(save_action)
        
        load_text = "Load" if lang == 'en' else "Laden" if lang == 'de' else "Завантажити" if lang == 'uk' else "Charger" if lang == 'fr' else "Load"
        load_action = QAction(load_text, self)
        load_action.triggered.connect(self.load_data)
        load_action.setShortcut("Ctrl+O")
        file_menu.addAction(load_action)

        close_all_text = "Close All Tabs" if lang == 'en' else "Alle Tabs schließen" if lang == 'de' else "Закрити всі вкладки" if lang == 'uk' else "Fermer tous les onglets" if lang == 'fr' else "Close All Tabs"
        close_all_action = QAction(close_all_text, self)
        close_all_action.triggered.connect(self.close_all_tabs)
        file_menu.addAction(close_all_action)

        file_menu.addSeparator()
        
        exit_text = "Exit" if lang == 'en' else "Beenden" if lang == 'de' else "Вихід" if lang == 'uk' else "Quitter" if lang == 'fr' else "Exit"
        exit_action = QAction(exit_text, self)
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)
        
        # меню налаштувань
        options_menu_text = "Settings" if lang == 'en' else "Einstellungen" if lang == 'de' else "Налаштування" if lang == 'uk' else "Paramètres" if lang == 'fr' else "Settings"
        options_menu = menubar.addMenu(options_menu_text)
        
        # створення підменю
        viz_menu_text = "Visualization Styles" if lang == 'en' else "Visualisierungsstile" if lang == 'de' else "Стилі візуалізації" if lang == 'uk' else "Styles de visualisation" if lang == 'fr' else "Visualization Styles"
        viz_menu = options_menu.addMenu(viz_menu_text)
        
        stick_text = "Stick (stick model)" if lang == 'en' else "Stick (Stabmodell)" if lang == 'de' else "Stick (паличкова модель)" if lang == 'uk' else "Stick (modèle en bâton)" if lang == 'fr' else "Stick (stick model)"
        stick_action = QAction(stick_text, self)
        stick_action.triggered.connect(lambda: self.set_visualization_style('stick'))
        viz_menu.addAction(stick_action)
        
        sphere_text = "Sphere (ball model)" if lang == 'en' else "Sphere (Kugelmodell)" if lang == 'de' else "Sphere (кулькова модель)" if lang == 'uk' else "Sphere (modèle en boule)" if lang == 'fr' else "Sphere (ball model)"
        sphere_action = QAction(sphere_text, self)
        sphere_action.triggered.connect(lambda: self.set_visualization_style('sphere'))
        viz_menu.addAction(sphere_action)
        
        line_text = "Line (wire model)" if lang == 'en' else "Line (Drahtmodell)" if lang == 'de' else "Line (дротова модель)" if lang == 'uk' else "Line (modèle en fil)" if lang == 'fr' else "Line (wire model)"
        line_action = QAction(line_text, self)
        line_action.triggered.connect(lambda: self.set_visualization_style('line'))
        viz_menu.addAction(line_action)
        
        viz_menu.addSeparator()
        
        show_h_text = "Show hydrogen atoms" if lang == 'en' else "Wasserstoffatome anzeigen" if lang == 'de' else "Показувати атоми водню" if lang == 'uk' else "Afficher les atomes d'hydrogène" if lang == 'fr' else "Show hydrogen atoms"
        self.show_hydrogens_action = QAction(show_h_text, self)
        self.show_hydrogens_action.triggered.connect(self.toggle_hydrogens)
        viz_menu.addAction(self.show_hydrogens_action)

        theme_menu_text = "Themes" if lang == 'en' else "Themen" if lang == 'de' else "Теми" if lang == 'uk' else "Thèmes" if lang == 'fr' else "Themes"
        theme_menu = options_menu.addMenu(theme_menu_text)
        
        light_theme_text = "Light Theme" if lang == 'en' else "Helles Theme" if lang == 'de' else "Світла тема" if lang == 'uk' else "Thème clair" if lang == 'fr' else "Light Theme"
        light_theme_action = QAction(light_theme_text, self)
        light_theme_action.triggered.connect(lambda: self.set_light_theme(True))
        theme_menu.addAction(light_theme_action)
        
        dark_theme_text = "Dark Theme" if lang == 'en' else "Dunkles Theme" if lang == 'de' else "Темна тема" if lang == 'uk' else "Thème sombre" if lang == 'fr' else "Dark Theme"
        dark_theme_action = QAction(dark_theme_text, self)
        dark_theme_action.triggered.connect(lambda: self.set_dark_theme(True))
        theme_menu.addAction(dark_theme_action)
        
        gray_theme_text = "Gray Theme" if lang == 'en' else "Graues Theme" if lang == 'de' else "Сіра тема" if lang == 'uk' else "Thème gris" if lang == 'fr' else "Gray Theme"
        gray_theme_action = QAction(gray_theme_text, self)
        gray_theme_action.triggered.connect(lambda: self.set_gray_theme(True))
        theme_menu.addAction(gray_theme_action)
        
        tab_position_text = "Tab Position" if lang == 'en' else "Tab-Position" if lang == 'de' else "Позиція вкладок" if lang == 'uk' else "Position des onglets" if lang == 'fr' else "Tab Position"
        tab_position_menu = options_menu.addMenu(tab_position_text)
        
        tab_north_text = "Tabs on Top" if lang == 'en' else "Tabs oben" if lang == 'de' else "Вкладки зверху" if lang == 'uk' else "Onglets en haut" if lang == 'fr' else "Tabs on Top"
        tab_north_action = QAction(tab_north_text, self)
        tab_north_action.triggered.connect(lambda: self.set_tab_position(QTabWidget.North))
        tab_position_menu.addAction(tab_north_action)
        
        tab_south_text = "Tabs on Bottom" if lang == 'en' else "Tabs unten" if lang == 'de' else "Вкладки знизу" if lang == 'uk' else "Onglets en bas" if lang == 'fr' else "Tabs on Bottom"
        tab_south_action = QAction(tab_south_text, self)
        tab_south_action.triggered.connect(lambda: self.set_tab_position(QTabWidget.South))
        tab_position_menu.addAction(tab_south_action)
        
        tab_west_text = "Tabs on Left" if lang == 'en' else "Tabs links" if lang == 'de' else "Вкладки зліва" if lang == 'uk' else "Onglets à gauche" if lang == 'fr' else "Tabs on Left"
        tab_west_action = QAction(tab_west_text, self)
        tab_west_action.triggered.connect(lambda: self.set_tab_position(QTabWidget.West))
        tab_position_menu.addAction(tab_west_action)
        
        tab_east_text = "Tabs on Right" if lang == 'en' else "Tabs rechts" if lang == 'de' else "Вкладки справа" if lang == 'uk' else "Onglets à droite" if lang == 'fr' else "Tabs on Right"
        tab_east_action = QAction(tab_east_text, self)
        tab_east_action.triggered.connect(lambda: self.set_tab_position(QTabWidget.East))
        tab_position_menu.addAction(tab_east_action)
        
        language_menu_text = "Language" if lang == 'en' else "Sprache" if lang == 'de' else "Мова" if lang == 'uk' else "Langue" if lang == 'fr' else "Language"
        language_menu = options_menu.addMenu(language_menu_text)
        
        english_action = QAction("English", self)
        english_action.triggered.connect(lambda: self.set_language("en"))
        language_menu.addAction(english_action)
        
        german_action = QAction("Deutsch", self)
        german_action.triggered.connect(lambda: self.set_language("de"))
        language_menu.addAction(german_action)
        
        ukrainian_action = QAction("Українська", self)
        ukrainian_action.triggered.connect(lambda: self.set_language("uk"))
        language_menu.addAction(ukrainian_action)
        
        french_action = QAction("Français", self)
        french_action.triggered.connect(lambda: self.set_language("fr"))
        language_menu.addAction(french_action)
        
        # допомога
        help_menu_text = "Help" if lang == 'en' else "Hilfe" if lang == 'de' else "Допомога" if lang == 'uk' else "Aide" if lang == 'fr' else "Help"
        help_menu = menubar.addMenu(help_menu_text)
        
        about_text = "About the Program" if lang == 'en' else "Über das Programm" if lang == 'de' else "Про програму" if lang == 'uk' else "À propos du programme" if lang == 'fr' else "About the Program"
        about_action = QAction(about_text, self)
        about_action.triggered.connect(self.show_about_dialog)
        help_menu.addAction(about_action)
        
        smiles_help_text = "How to Convert Formula to SMILES" if lang == 'en' else "Wie Formel in SMILES umwandeln" if lang == 'de' else "Як перетворити формулу в SMILES" if lang == 'uk' else "Comment convertir une formule en SMILES" if lang == 'fr' else "How to Convert Formula to SMILES"
        smiles_help_action = QAction(smiles_help_text, self)
        smiles_help_action.triggered.connect(self.show_smiles_help_dialog)
        help_menu.addAction(smiles_help_action)

    def set_tab_position(self, position):
        """Sets tab position"""
        self.tabs.setTabPosition(position)
    #загрузка тем
    def set_light_theme(self, save_setting=False):
        """Sets light theme"""
        QApplication.setStyle("Fusion")
        light_stylesheet = """
        QMainWindow, QDialog, QWidget {
            background-color: #f8f9fa;
            color: #212529;
        }
        QTabWidget::pane {
            border: 1px solid #dee2e6;
            background-color: #ffffff;
        }
        QTabBar::tab {
            background-color: #e9ecef;
            color: #495057;
            padding: 8px 16px;
            margin-right: 2px;
            border: 1px solid #dee2e6;
            border-bottom: none;
            border-top-left-radius: 4px;
            border-top-right-radius: 4px;
        }
        QTabBar::tab:selected {
            background-color: #ffffff;
            color: #495057;
            border-color: #dee2e6;
            border-bottom-color: #ffffff;
        }
        QTabBar::tab:hover:!selected {
            background-color: #f8f9fa;
        }
        QLineEdit, QTextEdit {
            background-color: #ffffff;
            color: #212529;
            border: 1px solid #ced4da;
            border-radius: 4px;
            padding: 6px 12px;
            font-size: 14px;
        }
        QLineEdit:focus, QTextEdit:focus {
            border-color: #3498db;
            outline: none;
        }
        QPushButton {
            background-color: #3498db;
            color: #ffffff;
            border: none;
            padding: 8px 16px;
            border-radius: 4px;
            font-weight: bold;
            font-size: 14px;
        }
        QPushButton:hover {
            background-color: #2980b9;
        }
        QPushButton:pressed {
            background-color: #21618c;
        }
        QTableWidget {
            background-color: #ffffff;
            color: #212529;
            gridline-color: #dee2e6;
            border: 1px solid #dee2e6;
            border-radius: 4px;
        }
        QHeaderView::section {
            background-color: #e9ecef;
            color: #495057;
            padding: 8px;
            border: none;
            font-weight: bold;
        }
        QMenuBar {
            background-color: #343a40;
            color: #ffffff;
        }
        QMenuBar::item {
            background-color: transparent;
            color: #ffffff;
            padding: 8px 12px;
        }
        QMenuBar::item:selected {
            background-color: #495057;
        }
        QMenu {
            background-color: #ffffff;
            color: #212529;
            border: 1px solid #dee2e6;
            border-radius: 4px;
        }
        QMenu::item {
            padding: 8px 24px 8px 12px;
        }
        QMenu::item:selected {
            background-color: #3498db;
            color: #ffffff;
        }
        QRadioButton {
            color: #212529;
            spacing: 8px;
        }
        QRadioButton::indicator {
            width: 16px;
            height: 16px;
            border-radius: 8px;
            border: 2px solid #6c757d;
        }
        QRadioButton::indicator:checked {
            background-color: #3498db;
            border-color: #3498db;
        }
        QLabel {
            color: #212529;
        }
        QScrollArea {
            background-color: #ffffff;
            border: none;
        }
        QStatusBar {
            background-color: #e9ecef;
            color: #495057;
        }
        """
        self.setStyleSheet(light_stylesheet)
        if save_setting:
            self.save_theme_settings("light")

    def set_dark_theme(self, save_setting=False):
        """Sets dark theme"""
        QApplication.setStyle("Fusion")
        dark_stylesheet = """
        QMainWindow, QDialog, QWidget {
            background-color: #2b2b2b;
            color: #ffffff;
        }
        QTabWidget::pane {
            border: 1px solid #555555;
            background-color: #3c3c3c;
        }
        QTabBar::tab {
            background-color: #404040;
            color: #cccccc;
            padding: 8px 16px;
            margin-right: 2px;
            border: 1px solid #555555;
            border-bottom: none;
            border-top-left-radius: 4px;
            border-top-right-radius: 4px;
        }
        QTabBar::tab:selected {
            background-color: #2b2b2b;
            color: #ffffff;
            border-color: #555555;
            border-bottom-color: #2b2b2b;
        }
        QTabBar::tab:hover:!selected {
            background-color: #353535;
        }
        QLineEdit, QTextEdit {
            background-color: #3c3c3c;
            color: #ffffff;
            border: 1px solid #555555;
            border-radius: 4px;
            padding: 6px 12px;
            font-size: 14px;
        }
        QLineEdit:focus, QTextEdit:focus {
            border-color: #3498db;
            outline: none;
        }
        QPushButton {
            background-color: #3498db;
            color: white;
            border: none;
            padding: 8px 16px;
            border-radius: 4px;
            font-weight: bold;
            font-size: 14px;
        }
        QPushButton:hover {
            background-color: #2980b9;
        }
        QPushButton:pressed {
            background-color: #21618c;
        }
        QTableWidget {
            background-color: #3c3c3c;
            color: #ffffff;
            gridline-color: #555555;
            border: 1px solid #555555;
            border-radius: 4px;
        }
        QHeaderView::section {
            background-color: #404040;
            color: #cccccc;
            padding: 8px;
            border: none;
            font-weight: bold;
        }
        QMenuBar {
            background-color: #343a40;
            color: #ffffff;
        }
        QMenuBar::item {
            background-color: transparent;
            color: #ffffff;
            padding: 8px 12px;
        }
        QMenuBar::item:selected {
            background-color: #495057;
        }
        QMenu {
            background-color: #3c3c3c;
            color: #ffffff;
            border: 1px solid #555555;
            border-radius: 4px;
        }
        QMenu::item {
            padding: 8px 24px 8px 12px;
        }
        QMenu::item:selected {
            background-color: #3498db;
            color: #ffffff;
        }
        QRadioButton {
            color: #ffffff;
            spacing: 8px;
        }
        QRadioButton::indicator {
            width: 16px;
            height: 16px;
            border-radius: 8px;
            border: 2px solid #cccccc;
        }
        QRadioButton::indicator:checked {
            background-color: #3498db;
            border-color: #3498db;
        }
        QLabel {
            color: #ffffff;
        }
        QScrollArea {
            background-color: #3c3c3c;
            border: none;
        }
        QStatusBar {
            background-color: #404040;
            color: #cccccc;
        }
        """
        self.setStyleSheet(dark_stylesheet)
        if save_setting:
            self.save_theme_settings("dark")

    def set_gray_theme(self, save_setting=False):
        """Sets gray theme"""
        QApplication.setStyle("Fusion")
        gray_stylesheet = """
        QMainWindow, QDialog, QWidget {
            background-color: #6c757d;
            color: #ffffff;
        }
        QTabWidget::pane {
            border: 1px solid #495057;
            background-color: #5a6268;
        }
        QTabBar::tab {
            background-color: #545b62;
            color: #ffffff;
            padding: 8px 16px;
            margin-right: 2px;
            border: 1px solid #495057;
            border-bottom: none;
            border-top-left-radius: 4px;
            border-top-right-radius: 4px;
        }
        QTabBar::tab:selected {
            background-color: #6c757d;
            color: #ffffff;
            border-color: #495057;
            border-bottom-color: #6c757d;
        }
        QTabBar::tab:hover:!selected {
            background-color: #5a6268;
        }
        QLineEdit, QTextEdit {
            background-color: #5a6268;
            color: #ffffff;
            border: 1px solid #495057;
            border-radius: 4px;
            padding: 6px 12px;
            font-size: 14px;
        }
        QLineEdit:focus, QTextEdit:focus {
            border-color: #3498db;
            outline: none;
        }
        QPushButton {
            background-color: #3498db;
            color: white;
            border: none;
            padding: 8px 16px;
            border-radius: 4px;
            font-weight: bold;
            font-size: 14px;
        }
        QPushButton:hover {
            background-color: #2980b9;
        }
        QPushButton:pressed {
            background-color: #21618c;
        }
        QTableWidget {
            background-color: #5a6268;
            color: #ffffff;
            gridline-color: #495057;
            border: 1px solid #495057;
            border-radius: 4px;
        }
        QHeaderView::section {
            background-color: #545b62;
            color: #ffffff;
            padding: 8px;
            border: none;
            font-weight: bold;
        }
        QMenuBar {
            background-color: #495057;
            color: #ffffff;
        }
        QMenuBar::item {
            background-color: transparent;
            color: #ffffff;
            padding: 8px 12px;
        }
        QMenuBar::item:selected {
            background-color: #5a6268;
        }
        QMenu {
            background-color: #5a6268;
            color: #ffffff;
            border: 1px solid #495057;
            border-radius: 4px;
        }
        QMenu::item {
            padding: 8px 24px 8px 12px;
        }
        QMenu::item:selected {
            background-color: #3498db;
            color: #ffffff;
        }
        QRadioButton {
            color: #ffffff;
            spacing: 8px;
        }
        QRadioButton::indicator {
            width: 16px;
            height: 16px;
            border-radius: 8px;
            border: 2px solid #cccccc;
        }
        QRadioButton::indicator:checked {
            background-color: #3498db;
            border-color: #3498db;
        }
        QLabel {
            color: #ffffff;
        }
        QScrollArea {
            background-color: #5a6268;
            border: none;
        }
        QStatusBar {
            background-color: #545b62;
            color: #ffffff;
        }
        """
        self.setStyleSheet(gray_stylesheet)
        if save_setting:
            self.save_theme_settings("gray")
    # допоміжні вікна
    def show_about_dialog(self):
        about_dialog = AboutDialog(self)
        about_dialog.exec_()
    
    def show_smiles_help_dialog(self):
        smiles_help_dialog = SmilesHelpDialog(self)
        smiles_help_dialog.exec_()

    def set_visualization_style(self, style):
        current_tab = self.tabs.currentWidget()
        if isinstance(current_tab, MainTabWidget):
            current_tab.set_visualization_style(style)

    def toggle_hydrogens(self):
        current_tab = self.tabs.currentWidget()
        if isinstance(current_tab, MainTabWidget):
            current_tab.toggle_hydrogens()
    # мова
    def set_language(self, lang):
        self.settings.setValue("language", lang)
        lang_current = self.current_language
        if lang_current == 'en':
            msg_title = "Language Change"
            msg_text = "Please restart the application to apply the language change."
        elif lang_current == 'de':
            msg_title = "Sprachänderung"
            msg_text = "Bitte starten Sie die Anwendung neu, um die Sprachänderung anzuwenden."
        elif lang_current == 'uk':
            msg_title = "Зміна мови"
            msg_text = "Будь ласка, перезапустіть додаток, щоб застосувати зміну мови."
        elif lang_current == 'fr':
            msg_title = "Changement de langue"
            msg_text = "Veuillez redémarrer l'application pour appliquer le changement de langue."
        QMessageBox.information(self, msg_title, msg_text)
    # збереження файлу
    def save_data(self):
        current_tab = self.tabs.currentWidget()
        if isinstance(current_tab, MainTabWidget):
            data = current_tab.get_data()
            data["current_tab"] = self.tabs.currentIndex()
            
            lang = self.current_language
            if lang == 'en':
                title = "Save File"
            elif lang == 'de':
                title = "Datei speichern"
            elif lang == 'uk':
                title = "Зберегти файл"
            elif lang == 'fr':
                title = "Enregistrer le fichier"
            file_name, _ = QFileDialog.getSaveFileName(self, title, "", "JSON Files (*.json);;All Files (*)")
            if file_name:
                with open(file_name, 'w') as f:
                    json.dump(data, f)
                if lang == 'en':
                    success_title = "Saved"
                    success_msg = "Data successfully saved!"
                elif lang == 'de':
                    success_title = "Gespeichert"
                    success_msg = "Daten erfolgreich gespeichert!"
                elif lang == 'uk':
                    success_title = "Збережено"
                    success_msg = "Дані успішно збережені!"
                elif lang == 'fr':
                    success_title = "Enregistré"
                    success_msg = "Données enregistrées avec succès !"
                QMessageBox.information(self, success_title, success_msg)
    # завантаження файлу
    def load_data(self):
        lang = self.current_language
        if lang == 'en':
            title = "Load File"
        elif lang == 'de':
            title = "Datei laden"
        elif lang == 'uk':
            title = "Завантажити файл"
        elif lang == 'fr':
            title = "Charger le fichier"
        file_name, _ = QFileDialog.getOpenFileName(self, title, "", "JSON Files (*.json);;All Files (*)")
        if file_name:
            with open(file_name, 'r') as f:
                data = json.load(f)
                
            current_tab = self.tabs.currentWidget()
            if isinstance(current_tab, MainTabWidget):
                current_tab.load_data(data)
            self.tabs.setCurrentIndex(data.get("current_tab", 0))
            if lang == 'en':
                success_title = "Loaded"
                success_msg = "Data successfully loaded!"
            elif lang == 'de':
                success_title = "Geladen"
                success_msg = "Daten erfolgreich geladen!"
            elif lang == 'uk':
                success_title = "Завантажено"
                success_msg = "Дані успішно завантажені!"
            elif lang == 'fr':
                success_title = "Chargé"
                success_msg = "Données chargées avec succès !"
            QMessageBox.information(self, success_title, success_msg)
    # експорт молекули
    def export_molecule_advanced(self):
        current_tab = self.tabs.currentWidget()
        if isinstance(current_tab, MainTabWidget):
            if not any([current_tab.current_mol, current_tab.current_mol_2, current_tab.current_mol_3]):
                return
                
            formats = "MOL Files (*.mol);;SDF Files (*.sdf);;SMILES Files (*.smi);;All Files (*)"
            lang = self.current_language
            if lang == 'en':
                title = "Export molecules"
            elif lang == 'de':
                title = "Moleküle exportieren"
            elif lang == 'uk':
                title = "Експорт молекул"
            elif lang == 'fr':
                title = "Exporter les molécules"
            filename, selected_filter = QFileDialog.getSaveFileName(
                self, title, "", formats)
                
            if filename:
                if selected_filter == "MOL Files (*.mol)":
                    self.export_as_mol(filename, current_tab)
                elif selected_filter == "SDF Files (*.sdf)":
                    self.export_as_sdf(filename, current_tab)
                elif selected_filter == "SMILES Files (*.smi)":
                    self.export_as_smiles(filename, current_tab)

    def export_as_mol(self, filename, tab):
        with open(filename, 'w') as f:
            molecules = [tab.current_mol, tab.current_mol_2, tab.current_mol_3]
            for mol in molecules:
                if mol:
                    mol_block = Chem.MolToMolBlock(mol)
                    f.write(mol_block)
                    f.write("\n$$$$\n")

    def export_as_sdf(self, filename, tab):
        writer = Chem.SDWriter(filename)
        molecules = [tab.current_mol, tab.current_mol_2, tab.current_mol_3]
        for mol in molecules:
            if mol:
                writer.write(mol)
        writer.close()

    def export_as_smiles(self, filename, tab):
        with open(filename, 'w') as f:
            molecules = [tab.current_mol, tab.current_mol_2, tab.current_mol_3]
            for mol in molecules:
                if mol:
                    smiles = Chem.MolToSmiles(mol)
                    f.write(smiles + "\n")

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event):
        files = [u.toLocalFile() for u in event.mimeData().urls()]
        for f in files:
            if f.endswith('.mol'):
                mol = Chem.MolFromMolFile(f)
                if mol:
                    smiles = Chem.MolToSmiles(mol)
                    current_tab = self.tabs.currentWidget()
                    if isinstance(current_tab, MainTabWidget):
                        current_tab.line_edit.setText(smiles)
                        current_tab.check_smiles()
            elif f.endswith('.sdf'):
                supplier = Chem.SDMolSupplier(f)
                for mol in supplier:
                    if mol:
                        smiles = Chem.MolToSmiles(mol)
                        current_tab = self.tabs.currentWidget()
                        if isinstance(current_tab, MainTabWidget):
                            current_tab.line_edit.setText(smiles)
                            current_tab.check_smiles()
                            break  # Take first molecule from SDF

if __name__ == "__main__":
    app = QApplication(sys.argv)
    

    app.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    app.setAttribute(Qt.AA_UseHighDpiPixmaps, True)
    
   
    app.setAttribute(Qt.AA_UseOpenGLES)
    

    icon_path = "icon.png"
    if os.path.exists(icon_path):
        app.setWindowIcon(QIcon(icon_path))
    
    window = MoleculeViewer()
    window.show()
    sys.exit(app.exec())
