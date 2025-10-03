# біблиотеки 😄😄😄
from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QLabel, QScrollArea, QFrame, QPushButton, QMessageBox, QApplication, QWidget, QFileDialog)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QIcon

class AboutDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(self.get_window_title())
        self.setFixedSize(900, 700)
        self.setWindowIcon(QIcon("icon.png"))
        
        
        self.is_dark_theme = self.is_dark_theme_active()
        
        self.apply_theme_styles()
        
        # головний віджет
        main_layout = QVBoxLayout()
        
        # назва
        title_label = QLabel(self.get_title_text())
        title_label.setStyleSheet(self.get_title_style())
        title_label.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(title_label)
        
        # прокрутка
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QFrame.NoFrame)
        scroll_area.setStyleSheet(self.get_scroll_area_style())
        
        content_widget = QLabel()
        content_layout = QVBoxLayout(content_widget)
        content_widget.setStyleSheet(self.get_content_style())
        
        # інформація
        info_text = self.get_info_text()
        
        info_label = QLabel(info_text)
        info_label.setWordWrap(True)
        info_label.setTextFormat(Qt.RichText)
        content_layout.addWidget(info_label)
        
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        separator.setStyleSheet(self.get_separator_style())
        content_layout.addWidget(separator)
        
        developer_info = self.get_developer_info()
        
        developer_label = QLabel(developer_info)
        developer_label.setAlignment(Qt.AlignRight)
        developer_label.setTextFormat(Qt.RichText)
        content_layout.addWidget(developer_label)
        
        scroll_area.setWidget(content_widget)
        main_layout.addWidget(scroll_area)
        
        # кнопка закритя
        close_button = QPushButton(self.get_close_button_text())
        close_button.setStyleSheet(self.get_button_style())
        close_button.clicked.connect(self.accept)
        main_layout.addWidget(close_button, alignment=Qt.AlignCenter)
        
        self.setLayout(main_layout)
    # вибір яка мова
    def get_window_title(self):
        lang = self.parent().current_language
        if lang == 'en':
            return "About the Program"
        elif lang == 'de':
            return "Über das Programm"
        elif lang == 'uk':
            return "Про програму"
        elif lang == 'fr':
            return "À propos du programme"
    
    def get_title_text(self):
        lang = self.parent().current_language
        if lang == 'en':
            return "ChemComparer - Tool for Molecule Analysis and Comparison"
        elif lang == 'de':
            return "ChemComparer - Werkzeug für Molekülanalyse und Vergleich"
        elif lang == 'uk':
            return "ChemComparer - Інструмент для аналізу та порівняння молекул"
        elif lang == 'fr':
            return "ChemComparer - Outil pour l'analyse et la comparaison des molécules"
    
    def get_close_button_text(self):
        lang = self.parent().current_language
        if lang == 'en':
            return "Close"
        elif lang == 'de':
            return "Schließen"
        elif lang == 'uk':
            return "Закрити"
        elif lang == 'fr':
            return "Fermer"
    
    def is_dark_theme_active(self):
        # проверка теми
        current_style = QApplication.style().objectName().lower()
        return 'dark' in current_style or 'fusion' in current_style and 'dark' in QApplication.palette().window().color().name().lower()
    
    def apply_theme_styles(self):
        """Applies styles according to the theme"""
        if self.is_dark_theme:
            self.setStyleSheet("""
                QDialog {
                    background-color: #2b2b2b;
                    color: #ffffff;
                }
                QLabel {
                    color: #ffffff;
                }
                QScrollArea {
                    background-color: #2b2b2b;
                    border: none;
                }
            """)
        else:
            self.setStyleSheet("""
                QDialog {
                    background-color: #ffffff;
                    color: #000000;
                }
                QLabel {
                    color: #000000;
                }
                QScrollArea {
                    background-color: #ffffff;
                    border: none;
                }
            """)
    
    def get_title_style(self):
        """Returns title style according to the theme"""
        if self.is_dark_theme:
            return "font-size: 20px; font-weight: bold; color: #3498db; padding: 10px;"
        else:
            return "font-size: 20px; font-weight: bold; color: #2c3e50; padding: 10px;"
    
    def get_scroll_area_style(self):
        """Returns scroll area style"""
        if self.is_dark_theme:
            return "background-color: #2b2b2b; border: none;"
        else:
            return "background-color: #ffffff; border: none;"
    
    def get_content_style(self):
        """Returns content style"""
        if self.is_dark_theme:
            return "background-color: #2b2b2b; color: #ffffff;"
        else:
            return "background-color: #ffffff; color: #000000;"
    
    def get_separator_style(self):
        """Returns separator style"""
        if self.is_dark_theme:
            return "background-color: #555555; margin: 15px 0;"
        else:
            return "background-color: #bdc3c7; margin: 15px 0;"
    
    def get_button_style(self):
        """Returns button style"""
        if self.is_dark_theme:
            return """
                QPushButton {
                    background-color: #3498db;
                    color: white;
                    border: none;
                    padding: 10px 20px;
                    font-weight: bold;
                    border-radius: 5px;
                }
                QPushButton:hover {
                    background-color: #2980b9;
                }
            """
        else:
            return """
                QPushButton {
                    background-color: #3498db;
                    color: white;
                    border: none;
                    padding: 10px 20px;
                    font-weight: bold;
                    border-radius: 5px;
                }
                QPushButton:hover {
                    background-color: #2980b9;
                }
            """
    # повернення теми і мови
    def get_info_text(self):
        text_color = "#e0e0e0" if self.is_dark_theme else "#333333"
        bg_color = "#3c3c3c" if self.is_dark_theme else "#f8f9fa"
        accent_color = "#3498db" if self.is_dark_theme else "#2c3e50"
        highlight_color = "#e74c3c" if self.is_dark_theme else "#e74c3c"
        
        lang = self.parent().current_language
        if lang == 'en':
            return f"""
            <div style='font-size: 12px; line-height: 1.6; color: {text_color};'>
            <p><strong style='color: {accent_color};'>ChemComparer Version 0.8.1</strong> - is a powerful tool that allows comprehensive
            analysis and comparison of molecular structures using SMILES notation.</p>
            
            <p style='color: {highlight_color}; font-weight: bold;'>🌟 This program is created to demonstrate the capabilities 
            of modern software in chemical informatics and molecular modeling.</p>
            
            <h3 style='color: {accent_color};'>🚀 New features in version 0.8.1:</h3>
            <ul>
                <li><strong>Tab system</strong> - work with multiple projects simultaneously</li>
                <li><strong>Theme support</strong> - light and dark interface themes</li>
                <li><strong>Drag & Drop</strong> - drag .mol and .sdf files directly into the program window</li>
                <li><strong>Hotkeys</strong> - speed up work with the program:
                    <ul>
                        <li>Ctrl+S - save project</li>
                        <li>Ctrl+O - load project</li>
                        <li>F11 - fullscreen mode</li>
                    </ul>
                </li>
                <li><strong>Improved error handling</strong> - more understandable messages</li>
                <li><strong>Automatic format detection</strong> - the program recognizes different file formats</li>
                <li><strong>Adaptive interface</strong> - adjusts to window size</li>
            </ul>
            
            <h3 style='color: {accent_color};'>📊 Main features:</h3>
            <ul>
                <li><strong>Molecule visualization</strong> in 2D and 3D formats with high quality</li>
                <li><strong>Comparative analysis</strong> of up to three molecules simultaneously</li>
                <li><strong>Advanced export</strong> in MOL, SDF, SMILES formats</li>
                <li><strong>Various display styles</strong> for 3D models (stick, ball, wire)</li>
                <li><strong>Flexible display management</strong> of hydrogen atoms</li>
                <li><strong>Interactive comparative diagrams</strong> for property visualization</li>
                <li><strong>Fullscreen mode</strong> for tables and 3D models</li>
            </ul>
            
            <h3 style='color: {accent_color};'>🎯 How to use the program:</h3>
            <ol>
                <li>Enter the molecule SMILES code in the main input field</li>
                <li>Click "Check SMILES and launch" to load the molecule</li>
                <li>For comparison, enter SMILES of other molecules in the corresponding fields</li>
                <li>Use the "Add second/third molecule" buttons</li>
                <li>Switch between 2D and 3D modes using radio buttons</li>
                <li>Use the "Settings" menu to change visualization style</li>
                <li>Export results in the desired format</li>
                <li>Create comparative diagrams for property analysis</li>
            </ol>
            
            <h3 style='color: {accent_color};'>📁 Supported formats:</h3>
            <ul>
                <li><strong>Input:</strong> SMILES notation, .mol, .sdf files</li>
                <li><strong>Export:</strong> MOL, SDF, SMILES</li>
                <li><strong>Session saving:</strong> JSON</li>
                <li><strong>Visualization:</strong> PNG (2D), interactive 3D</li>
            </ul>
            
            <h3 style='color: {accent_color};'>⚙️ Technology stack:</h3>
            <ul>
                <li><strong>Python 3</strong> - main programming language</li>
                <li><strong>PyQt5</strong> - graphical user interface</li>
                <li><strong>RDKit</strong> - chemical computations and analysis</li>
                <li><strong>py3Dmol</strong> - 3D molecule visualization</li>
                <li><strong>Matplotlib</strong> - creating diagrams and graphs</li>
                <li><strong>Pandas</strong> - data processing for diagrams</li>
            </ul>
            
            <div style='background-color: {bg_color}; padding: 15px; border-radius: 5px;'>
            <strong>🎓 Educational aspect:</strong> This project demonstrates the integration of modern software technologies 
            with chemical research. It includes work with graphical interface, chemical data processing, 
            3D visualization and creation of interactive diagrams.
            </div>
            
            <p style='color: #7f8c8d; font-style: italic;'>
            The program is developed as a personal educational project for myself, with code assisted by AI.
            All functions are implemented with emphasis on the learning process and demonstration of capabilities
            of modern software for chemists.
            </p>
            </div>
            """
        elif lang == 'de':
            return f"""
            <div style='font-size: 12px; line-height: 1.6; color: {text_color};'>
            <p><strong style='color: {accent_color};'>ChemComparer Version 0.8.1</strong> - ist ein leistungsstarkes Tool, das eine umfassende
            Analyse und Vergleich von Molekülstrukturen mit Hilfe der SMILES-Notation ermöglicht.</p>
            
            <p style='color: {highlight_color}; font-weight: bold;'>🌟 Dieses Programm wurde erstellt, um die Möglichkeiten 
            moderner Software in der chemischen Informatik und Molekülmodellierung zu demonstrieren.</p>
            
            <h3 style='color: {accent_color};'>🚀 Neue Funktionen in Version 0.8.1:</h3>
            <ul>
                <li><strong>Tab-System</strong> - arbeiten Sie gleichzeitig mit mehreren Projekten</li>
                <li><strong>Themenunterstützung</strong> - helle und dunkle Interface-Themen</li>
                <li><strong>Drag & Drop</strong> - ziehen Sie .mol- und .sdf-Dateien direkt in das Programmfenster</li>
                <li><strong>Tastenkürzel</strong> - beschleunigen Sie die Arbeit mit dem Programm:
                    <ul>
                        <li>Ctrl+S - Projekt speichern</li>
                        <li>Ctrl+O - Projekt laden</li>
                        <li>F11 - Vollbildmodus</li>
                    </ul>
                </li>
                <li><strong>Verbesserte Fehlerbehandlung</strong> - verständlichere Meldungen</li>
                <li><strong>Automatische Formaterkennung</strong> - das Programm erkennt verschiedene Dateiformate</li>
                <li><strong>Adaptives Interface</strong> - passt sich der Fenstergröße an</li>
            </ul>
            
            <h3 style='color: {accent_color};'>📊 Hauptfunktionen:</h3>
            <ul>
                <li><strong>Molekülvisualisierung</strong> in 2D- und 3D-Formaten mit hoher Qualität</li>
                <li><strong>Vergleichsanalyse</strong> von bis zu drei Molekülen gleichzeitig</li>
                <li><strong>Erweiterter Export</strong> in MOL-, SDF-, SMILES-Formaten</li>
                <li><strong>Verschiedene Anzeigestile</strong> für 3D-Modelle (Stab, Kugel, Draht)</li>
                <li><strong>Flexible Anzeigeverwaltung</strong> von Wasserstoffatomen</li>
                <li><strong>Interaktive Vergleichsdiagramme</strong> zur Visualisierung von Eigenschaften</li>
                <li><strong>Vollbildmodus</strong> für Tabellen und 3D-Modelle</li>
            </ul>
            
            <h3 style='color: {accent_color};'>🎯 So verwenden Sie das Programm:</h3>
            <ol>
                <li>Geben Sie den SMILES-Code des Moleküls in das Haupt-Eingabefeld ein</li>
                <li>Klicken Sie auf "SMILES überprüfen und starten", um das Molekül zu laden</li>
                <li>Für den Vergleich geben Sie SMILES anderer Moleküle in die entsprechenden Felder ein</li>
                <li>Verwenden Sie die Schaltflächen "Zweites/drittes Molekül hinzufügen"</li>
                <li>Wechseln Sie zwischen 2D- und 3D-Modi mit Hilfe der Radiobuttons</li>
                <li>Verwenden Sie das Menü "Einstellungen", um den Visualisierungsstil zu ändern</li>
                <li>Exportieren Sie die Ergebnisse im gewünschten Format</li>
                <li>Erstellen Sie Vergleichsdiagramme zur Analyse der Eigenschaften</li>
            </ol>
            
            <h3 style='color: {accent_color};'>📁 Unterstützte Formate:</h3>
            <ul>
                <li><strong>Eingabe:</strong> SMILES-Notation, .mol-, .sdf-Dateien</li>
                <li><strong>Export:</strong> MOL, SDF, SMILES</li>
                <li><strong>Sitzungsspeicherung:</strong> JSON</li>
                <li><strong>Visualisierung:</strong> PNG (2D), interaktives 3D</li>
            </ul>
            
            <h3 style='color: {accent_color};'>⚙️ Technologie-Stack:</h3>
            <ul>
                <li><strong>Python 3</strong> - Hauptsprache der Programmierung</li>
                <li><strong>PyQt5</strong> - grafische Benutzeroberfläche</li>
                <li><strong>RDKit</strong> - chemische Berechnungen und Analyse</li>
                <li><strong>py3Dmol</strong> - 3D-Molekülvisualisierung</li>
                <li><strong>Matplotlib</strong> - Erstellung von Diagrammen und Graphen</li>
                <li><strong>Pandas</strong> - Datenverarbeitung für Diagramme</li>
            </ul>
            
            <div style='background-color: {bg_color}; padding: 15px; border-radius: 5px;'>
            <strong>🎓 Bildungsaspekt:</strong> Dieses Projekt demonstriert die Integration moderner Softwaretechnologien 
            mit chemischer Forschung. Es umfasst Arbeit mit grafischer Oberfläche, Verarbeitung chemischer Daten, 
            3D-Visualisierung und Erstellung interaktiver Diagramme.
            </div>
            
            <p style='color: #7f8c8d; font-style: italic;'>
            Das Programm wurde als persönliches Bildungsprojekt für mich selbst entwickelt, mit Code-Unterstützung durch KI.
            Alle Funktionen sind mit Schwerpunkt auf dem Lernprozess und der Demonstration der Möglichkeiten
            moderner Software für Chemiker implementiert.
            </p>
            </div>
            """
        elif lang == 'uk':
            return f"""
            <div style='font-size: 12px; line-height: 1.6; color: {text_color};'>
            <p><strong style='color: {accent_color};'>ChemComparer Версія 0.8.1</strong> - це потужний інструмент, який дозволяє комплексний
            аналіз та порівняння молекулярних структур за допомогою SMILES нотації.</p>
            
            <p style='color: {highlight_color}; font-weight: bold;'>🌟 Ця програма створена для демонстрації можливостей 
            сучасного програмного забезпечення в хімічній інформатиці та молекулярному моделюванні.</p>
            
            <h3 style='color: {accent_color};'>🚀 Нові функції в версії 0.8.1:</h3>
            <ul>
                <li><strong>Система вкладок</strong> - працюйте з кількома проектами одночасно</li>
                <li><strong>Підтримка тем</strong> - світла та темна теми інтерфейсу</li>
                <li><strong>Drag & Drop</strong> - перетягуйте файли .mol та .sdf прямо у вікно програми</li>
                <li><strong>Гарячі клавіші</strong> - прискорьте роботу з програмою:
                    <ul>
                        <li>Ctrl+S - зберегти проект</li>
                        <li>Ctrl+O - завантажити проект</li>
                        <li>F11 - повноекранний режим</li>
                    </ul>
                </li>
                <li><strong>Покращена обробка помилок</strong> - більш зрозумілі повідомлення</li>
                <li><strong>Автоматичне визначення формату</strong> - програма розпізнає різні формати файлів</li>
                <li><strong>Адаптивний інтерфейс</strong> - підлаштовується під розмір вікна</li>
            </ul>
            
            <h3 style='color: {accent_color};'>📊 Основні можливості:</h3>
            <ul>
                <li><strong>Візуалізація молекул</strong> у 2D та 3D форматах з високою якістю</li>
                <li><strong>Порівняльний аналіз</strong> до трьох молекул одночасно</li>
                <li><strong>Розширений експорт</strong> у форматах MOL, SDF, SMILES</li>
                <li><strong>Різноманітні стилі відображення</strong> 3D моделей (паличковий, кульковий, дротовий)</li>
                <li><strong>Гнучке керування відображенням</strong> атомів водню</li>
                <li><strong>Інтерактивні порівняльні діаграми</strong> для візуалізації властивостей</li>
                <li><strong>Повноекранний режим</strong> для таблиць та 3D моделей</li>
            </ul>
            
            <h3 style='color: {accent_color};'>🎯 Як користуватися програмою:</h3>
            <ol>
                <li>Введіть SMILES-код молекули у головне поле введення</li>
                <li>Натисніть "Перевірити SMILES і запустити" для завантаження молекули</li>
                <li>Для порівняння введіть SMILES других молекул у відповідні поля</li>
                <li>Використовуйте кнопки "Додати другу/третю молекулу"</li>
                <li>Перемикайтеся між 2D та 3D режимами за допомогою радіокнопок</li>
                <li>Використовуйте меню "Налаштування" для зміни стилю візуалізації</li>
                <li>Експортуйте результати у потрібному форматі</li>
                <li>Створюйте порівняльні діаграми для аналізу властивостей</li>
            </ol>
            
            <h3 style='color: {accent_color};'>📁 Підтримувані формати:</h3>
            <ul>
                <li><strong>Введення:</strong> SMILES нотація, файли .mol, .sdf</li>
                <li><strong>Експорт:</strong> MOL, SDF, SMILES</li>
                <li><strong>Збереження сеансів:</strong> JSON</li>
                <li><strong>Візуалізація:</strong> PNG (2D), інтерактивний 3D</li>
            </ul>
            
            <h3 style='color: {accent_color};'>⚙️ Технологічний стек:</h3>
            <ul>
                <li><strong>Python 3</strong> - основна мова програмування</li>
                <li><strong>PyQt5</strong> - графічний інтерфейс користувача</li>
                <li><strong>RDKit</strong> - хімічні обчислення та аналіз</li>
                <li><strong>py3Dmol</strong> - 3D візуалізація молекул</li>
                <li><strong>Matplotlib</strong> - створення діаграм та графіків</li>
                <li><strong>Pandas</strong> - обробка даних для діаграм</li>
            </ul>
            
            <div style='background-color: {bg_color}; padding: 15px; border-radius: 5px;'>
            <strong>🎓 Навчальний аспект:</strong> Цей проект демонструє інтеграцію сучасних програмних технологій 
            з хімічними дослідженнями. Він включає роботу з графічним інтерфейсом, обробку хімічних даних, 
            3D візуалізацію та створення інтерактивних діаграм.
            </div>
            
            <p style='color: #7f8c8d; font-style: italic;'>
            Програма розроблена як особистий навчальний проект для себе, з кодом допомагав ШІ.
            Всі функції реалізовані з акцентом на навчальний процес та демонстрацію можливостей
            сучасного програмного забезпечення для хіміків.
            </p>
            </div>
            """
        elif lang == 'fr':
            return f"""
            <div style='font-size: 12px; line-height: 1.6; color: {text_color};'>
            <p><strong style='color: {accent_color};'>ChemComparer Version 0.8.1</strong> - est un outil puissant qui permet une analyse complète
            et comparaison des structures moléculaires à l'aide de la notation SMILES.</p>
            
            <p style='color: {highlight_color}; font-weight: bold;'>🌟 Ce programme est créé pour démontrer les capacités 
            du logiciel moderne en informatique chimique et modélisation moléculaire.</p>
            
            <h3 style='color: {accent_color};'>🚀 Nouvelles fonctionnalités dans la version 0.8.1 :</h3>
            <ul>
                <li><strong>Système d'onglets</strong> - travaillez avec plusieurs projets simultanément</li>
                <li><strong>Support des thèmes</strong> - thèmes clair et sombre de l'interface</li>
                <li><strong>Drag & Drop</strong> - glissez les fichiers .mol et .sdf directement dans la fenêtre du programme</li>
                <li><strong>Raccourcis clavier</strong> - accélerez le travail avec le programme :
                    <ul>
                        <li>Ctrl+S - sauvegarder le projet</li>
                        <li>Ctrl+O - charger le projet</li>
                        <li>F11 - mode plein écran</li>
                    </ul>
                </li>
                <li><strong>Gestion des erreurs améliorée</strong> - messages plus compréhensibles</li>
                <li><strong>Détection automatique du format</strong> - le programme reconnaît différents formats de fichiers</li>
                <li><strong>Interface adaptive</strong> - s'ajuste à la taille de la fenêtre</li>
            </ul>
            
            <h3 style='color: {accent_color};'>📊 Fonctionnalités principales :</h3>
            <ul>
                <li><strong>Visualisation des molécules</strong> aux formats 2D et 3D avec haute qualité</li>
                <li><strong>Analyse comparative</strong> jusqu'à trois molécules simultanément</li>
                <li><strong>Export avancé</strong> aux formats MOL, SDF, SMILES</li>
                <li><strong>Styles d'affichage variés</strong> pour les modèles 3D (bâton, boule, fil)</li>
                <li><strong>Gestion flexible de l'affichage</strong> des atomes d'hydrogène</li>
                <li><strong>Diagrammes comparatifs interactifs</strong> pour la visualisation des propriétés</li>
                <li><strong>Mode plein écran</strong> pour les tables et les modèles 3D</li>
            </ul>
            
            <h3 style='color: {accent_color};'>🎯 Comment utiliser le programme :</h3>
            <ol>
                <li>Entrez le code SMILES de la molécule dans le champ d'entrée principal</li>
                <li>Cliquez sur "Vérifier SMILES et lancer" pour charger la molécule</li>
                <li>Pour la comparaison, entrez SMILES d'autres molécules dans les champs correspondants</li>
                <li>Utilisez les boutons "Ajouter deuxième/troisième molécule"</li>
                <li>Basculer entre les modes 2D et 3D à l'aide des boutons radio</li>
                <li>Utilisez le menu "Paramètres" pour changer le style de visualisation</li>
                <li>Exportez les résultats au format souhaité</li>
                <li>Créez des diagrammes comparatifs pour l'analyse des propriétés</li>
            </ol>
            
            <h3 style='color: {accent_color};'>📁 Formats pris en charge :</h3>
            <ul>
                <li><strong>Entrée :</strong> Notation SMILES, fichiers .mol, .sdf</li>
                <li><strong>Export :</strong> MOL, SDF, SMILES</li>
                <li><strong>Sauvegarde de session :</strong> JSON</li>
                <li><strong>Visualisation :</strong> PNG (2D), 3D interactif</li>
            </ul>
            
            <h3 style='color: {accent_color};'>⚙️ Pile technologique :</h3>
            <ul>
                <li><strong>Python 3</strong> - langage de programmation principal</li>
                <li><strong>PyQt5</strong> - interface utilisateur graphique</li>
                <li><strong>RDKit</strong> - calculs chimiques et analyse</li>
                <li><strong>py3Dmol</strong> - visualisation 3D des molécules</li>
                <li><strong>Matplotlib</strong> - création de diagrammes et graphiques</li>
                <li><strong>Pandas</strong> - traitement des données pour les diagrammes</li>
            </ul>
            
            <div style='background-color: {bg_color}; padding: 15px; border-radius: 5px;'>
            <strong>🎓 Aspect éducatif :</strong> Ce projet démontre l'intégration des technologies logicielles modernes 
            avec la recherche chimique. Il inclut le travail avec l'interface graphique, le traitement des données chimiques, 
            la visualisation 3D et la création de diagrammes interactifs.
            </div>
            
            <p style='color: #7f8c8d; font-style: italic;'>
            Le programme est développé comme un projet éducatif personnel pour moi-même, avec le code assisté par IA.
            Toutes les fonctions sont implémentées avec un accent sur le processus d'apprentissage et la démonstration des capacités
            du logiciel moderne pour les chimistes.
            </p>
            </div>
            """

    def get_developer_info(self):
        text_color = "#cccccc" if self.is_dark_theme else "#7f8c8d"
        
        lang = self.parent().current_language
        if lang == 'en':
            return f"""
            <div style='font-size: 11px; color: {text_color};'>
            <p><strong>Developer:</strong> Daniil Tovkach</p>
            <p><strong>Version:</strong> 0.8.1</p>
            <p><strong>Release Date:</strong> Current</p>
            <p><strong>License:</strong> MIT Open Source</p>
            <p><strong>Contact:</strong> dtovkac4@gmail.com</p>
            <p><strong>GitHub:</strong> -</p>
            </div>
            """
        elif lang == 'de':
            return f"""
            <div style='font-size: 11px; color: {text_color};'>
            <p><strong>Entwickler:</strong> Daniil Tovkach</p>
            <p><strong>Version:</strong> 0.8.1</p>
            <p><strong>Veröffentlichungsdatum:</strong> Aktuell</p>
            <p><strong>Lizenz:</strong> MIT Open Source</p>
            <p><strong>Kontakt:</strong> dtovkac4@gmail.com</p>
            <p><strong>GitHub:</strong> -</p>
            </div>
            """
        elif lang == 'uk':
            return f"""
            <div style='font-size: 11px; color: {text_color};'>
            <p><strong>Розробник:</strong> Данііл Товкач</p>
            <p><strong>Версія:</strong> 0.8.1</p>
            <p><strong>Дата релізу:</strong> Поточна</p>
            <p><strong>Ліцензія:</strong> MIT Open Source</p>
            <p><strong>Зв'язок:</strong> dtovkac4@gmail.com</p>
            <p><strong>GitHub:</strong> -</p>
            </div>
            """
        elif lang == 'fr':
            return f"""
            <div style='font-size: 11px; color: {text_color};'>
            <p><strong>Développeur:</strong> Daniil Tovkach</p>
            <p><strong>Version:</strong> 0.8.1</p>
            <p><strong>Date de sortie:</strong> Actuelle</p>
            <p><strong>Licence:</strong> MIT Open Source</p>
            <p><strong>Contact:</strong> dtovkac4@gmail.com</p>
            <p><strong>GitHub:</strong> -</p>
            </div>
            """
# те саме створення але для SmilesHelpDialog і налаштування мови
class SmilesHelpDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(self.get_window_title())
        self.setFixedSize(950, 750)
        self.setWindowIcon(QIcon("icon.png"))
        

        main_layout = QVBoxLayout()
        

        title_label = QLabel(self.get_title_text())
        title_label.setStyleSheet("font-size: 20px; font-weight: bold; color: #2c3e50; padding: 10px;")
        title_label.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(title_label)
        

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QFrame.NoFrame)
        

        content_widget = QWidget()
        content_layout = QVBoxLayout(content_widget)
        

        help_text = self.get_help_text()
        
        help_label = QLabel(help_text)
        help_label.setWordWrap(True)
        help_label.setTextFormat(Qt.RichText)
        help_label.setOpenExternalLinks(True)  
        content_layout.addWidget(help_label)
        

        scroll_area.setWidget(content_widget)
        main_layout.addWidget(scroll_area)
        

        close_button = QPushButton(self.get_close_button_text())
        close_button.setStyleSheet("""
            QPushButton {
                background-color: #3498db;
                color: white;
                border: none;
                padding: 10px 20px;
                font-weight: bold;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #2980b9;
            }
        """)
        close_button.clicked.connect(self.accept)
        main_layout.addWidget(close_button, alignment=Qt.AlignCenter)
        
        self.setLayout(main_layout)

    def get_window_title(self):
        lang = self.parent().current_language
        if lang == 'en':
            return "SMILES Notation Help"
        elif lang == 'de':
            return "Hilfe zur SMILES-Notation"
        elif lang == 'uk':
            return "Довідка з нотації SMILES"
        elif lang == 'fr':
            return "Aide à la notation SMILES"
    
    def get_title_text(self):
        lang = self.parent().current_language
        if lang == 'en':
            return "SMILES - Simplified Molecular Input Line Entry System"
        elif lang == 'de':
            return "SMILES - Simplified Molecular Input Line Entry System"
        elif lang == 'uk':
            return "SMILES - Simplified Molecular Input Line Entry System"
        elif lang == 'fr':
            return "SMILES - Simplified Molecular Input Line Entry System"
    
    def get_close_button_text(self):
        lang = self.parent().current_language
        if lang == 'en':
            return "Close"
        elif lang == 'de':
            return "Schließen"
        elif lang == 'uk':
            return "Закрити"
        elif lang == 'fr':
            return "Fermer"
    
    def get_help_text(self):
        lang = self.parent().current_language
        if lang == 'en':
            return """
            <div style='font-size: 12px; line-height: 1.6;'>
            <h3 style='color: #3498db;'>What is SMILES?</h3>
            <p><strong>SMILES</strong> - is a standard text format for representing chemical structures, 
            developed for convenient input and exchange of chemical information.</p>
            
            <h3 style='color: #3498db;'>Basic syntax rules:</h3>
            <ol>
                <li><strong>Atoms:</strong> Designated by chemical symbols (C, O, N, Cl, Na etc.)
                    <ul>
                        <li>Uppercase letters: C, O, N - organic atoms</li>
                        <li>Lowercase letters: c, o, n - aromatic atoms</li>
                        <li>Elements in brackets: [Na+], [Cl-], [Fe++] - ions and special atoms</li>
                    </ul>
                </li>
                
                <li><strong>Bonds:</strong>
                    <ul>
                        <li>Single bond: not indicated or '-' (C-C or CC)</li>
                        <li>Double bond: '=' (C=O - carbonyl group)</li>
                        <li>Triple bond: '#' (C#N - nitrile group)</li>
                        <li>Aromatic bonds: automatically determined for lowercase letters</li>
                    </ul>
                </li>
                
                <li><strong>Rings:</strong> Numbers are used to denote connection points
                    <ul>
                        <li>C1CCCC1 - cyclopentane</li>
                        <li>c1ccccc1 - benzene (aromatic)</li>
                        <li>N1CCOCC1 - morpholine</li>
                    </ul>
                </li>
                
                <li><strong>Branching:</strong> Denoted by parentheses
                    <ul>
                        <li>CC(=O)O - acetic acid</li>
                        <li>CC(C)C - isobutane</li>
                        <li>NC(=O)C1CC1 - cyclopropanecarboxamide</li>
                    </ul>
                </li>
                
                <li><strong>Stereochemistry:</strong> Denoted by special symbols
                    <ul>
                        <li>C/C=C/C - trans-butene-2</li>
                        <li>C/C=C\\C - cis-butene-2</li>
                        <li>C[@H](F)(Cl)Br - chiral center</li>
                    </ul>
                </li>
            </ol>
            
            <h3 style='color: #3498db;'>Detailed conversion examples:</h3>
            
            <table border='1' cellpadding='8' style='border-collapse: collapse; width: 100%; margin: 15px 0;'>
            <tr style='background-color: #f8f9fa;'>
                <th style='padding: 10px;'>Molecule</th>
                <th style='padding: 10px;'>SMILES</th>
                <th style='padding: 10px;'>Explanation</th>
            </tr>
            <tr>
                <td style='padding: 8px;'><strong>Water</strong></td>
                <td style='padding: 8px;'><code>O</code></td>
                <td style='padding: 8px;'>Simple molecule - just oxygen symbol</td>
            </tr>
            <tr>
                <td style='padding: 8px;'><strong>Ethanol</strong></td>
                <td style='padding: 8px;'><code>CCO</code></td>
                <td style='padding: 8px;'>Chain of 2 carbons and oxygen</td>
            </tr>
            <tr>
                <td style='padding: 8px;'><strong>Benzene</strong></td>
                <td style='padding: 8px;'><code>c1ccccc1</code></td>
                <td style='padding: 8px;'>Aromatic ring (lowercase letters)</td>
            </tr>
            <tr>
                <td style='padding: 8px;'><strong>Aspirin</strong></td>
                <td style='padding: 8px;'><code>CC(=O)OC1=CC=CC=C1C(=O)O</code></td>
                <td style='padding: 8px;'>Complex structure with functional groups</td>
            </tr>
            <tr>
                <td style='padding: 8px;'><strong>Caffeine</strong></td>
                <td style='padding: 8px;'><code>CN1C=NC2=C1C(=O)N(C(=O)N2C)C</code></td>
                <td style='padding: 8px;'>Heterocyclic rings with methyl groups</td>
            </tr>
            <tr>
                <td style='padding: 8px;'><strong>Glucose</strong></td>
                <td style='padding: 8px;'><code>C(C1C(C(C(C(O1)O)O)O)O)O</code></td>
                <td style='padding: 8px;'>Complex structure with many OH groups</td>
            </tr>
            </table>
            
            <h3 style='color: #3498db;'>Complex example parsing:</h3>
            <p><strong>Molecule:</strong> Acetylsalicylic acid (Aspirin)</p>
            <p><strong>SMILES:</strong> <code style='background-color: #f8f9fa; padding: 5px;'>CC(=O)OC1=CC=CC=C1C(=O)O</code></p>
            
            <p><strong>Parsing by parts:</strong></p>
            <ul>
                <li><code>CC</code> - methyl group (CH3-)</li>
                <li><code>(=O)</code> - carbonyl group (C=O)</li>
                <li><code>O</code> - oxygen of ester bond</li>
                <li><code>c1ccccc1</code> - benzene ring (aromatic)</li>
                <li><code>C(=O)O</code> - carboxyl group (-COOH)</li>
            </ul>
            
            <h3 style='color: #3498db;'>Special symbols and their meanings:</h3>
            <ul>
                <li><strong>[]</strong> - atoms with specific properties: 
                    <a href="https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html">[Na+]</a>, 
                    <a href="https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system#Atoms">[Fe++]</a>, 
                    <a href="https://rdkit.org/docs/GettingStartedInPython.html">[nH]</a>
                </li>
                <li><strong>@</strong> - stereochemistry: 
                    <a href="https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html">@</a>, 
                    <a href="https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system#Stereochemistry">@@</a>
                </li>
                <li><strong>.</strong> - molecule separation: 
                    <a href="https://rdkit.org/docs/Cookbook.html">Na.Cl</a> (salt)
                </li>
                <li><strong>/</strong> and <strong>\\</strong> - stereochemical bonds: 
                    <a href="https://www.cheminfo.org/flavor/malaria/Utilities/SMILES_canonization/index.html">C/C=C/C</a>
                </li>
            </ul>
            
            <h3 style='color: #3498db;'>Useful resources and tools:</h3>
            <ul>
                <li><a href="https://pubchem.ncbi.nlm.nih.gov/">PubChem</a> - large database of chemical compounds</li>
                <li><a href="https://www.chemspider.com/">ChemSpider</a> - search by chemical structures</li>
                <li><a href="https://www.rdkit.org/">RDKit</a> - library for chemical informatics</li>
                <li><a href="https://cactus.nci.nih.gov/chemical/structure">NCI Cactus</a> - converter of names to SMILES</li>
                <li><a href="https://www.cheminfo.org/flavor/malaria/Utilities/SMILES_canonization/index.html">SMILES Canonizer</a> - validation and canonization</li>
                <li><a href="https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system">Wikipedia SMILES</a> - detailed documentation</li>
                <li><a href="https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html">Daylight Theory</a> - official specification</li>
            </ul>
            
            <h3 style='color: #3498db;'>How to get SMILES for your molecules:</h3>
            <ul>
                <li><strong>Chemical editors:</strong> 
                    <a href="https://chemaxon.com/products/marvin">MarvinSketch</a>, 
                    <a href="https://www.perkinelmer.com/category/chemdraw">ChemDraw</a>
                </li>
                <li><strong>Online tools:</strong> 
                    <a href="https://www.molview.org/">MolView</a>, 
                    <a href="https://www.cheminfo.org/">ChemDoodle</a>
                </li>
                <li><strong>Databases:</strong> 
                    <a href="https://www.drugbank.ca/">DrugBank</a>, 
                    <a href="https://www.ebi.ac.uk/chembl/">ChEMBL</a>
                </li>
                <li><strong>Mobile apps:</strong> 
                    <a href="https://play.google.com/store/apps/details?id=com.cheminfomatics.molprime">MolPrime+</a>
                </li>
            </ul>
            
            <div style='background-color: #e8f5e8; padding: 15px; border-radius: 5px; margin: 15px 0;'>
            <strong>💡 Tip:</strong> For complex molecules, we recommend using graphical editors, 
            which automatically generate SMILES code from the drawn structure.
            </div>
            
            <p style='color: #7f8c8d; font-style: italic;'>
            This help is created to assist in learning SMILES notation. For deeper study, 
            we recommend referring to official documentation and specialized resources.
            </p>
            </div>
            """
        elif lang == 'de':
            return """
            <div style='font-size: 12px; line-height: 1.6;'>
            <h3 style='color: #3498db;'>Was ist SMILES?</h3>
            <p><strong>SMILES</strong> - ist ein standardisiertes Textformat zur Darstellung chemischer Strukturen, 
            entwickelt für bequeme Eingabe und Austausch chemischer Informationen.</p>
            
            <h3 style='color: #3498db;'>Grundlegende Syntaxregeln:</h3>
            <ol>
                <li><strong>Atome:</strong> Werden durch chemische Symbole bezeichnet (C, O, N, Cl, Na usw.)
                    <ul>
                        <li>Großbuchstaben: C, O, N - organische Atome</li>
                        <li>Kleinbuchstaben: c, o, n - aromatische Atome</li>
                        <li>Elemente in Klammern: [Na+], [Cl-], [Fe++] - Ionen und spezielle Atome</li>
                    </ul>
                </li>
                
                <li><strong>Bindungen:</strong>
                    <ul>
                        <li>Einfache Bindung: nicht angegeben oder '-' (C-C oder CC)</li>
                        <li>Doppelbindung: '=' (C=O - Carbonylgruppe)</li>
                        <li>Dreifachbindung: '#' (C#N - Nitrilgruppe)</li>
                        <li>Aromatische Bindungen: automatisch für Kleinbuchstaben bestimmt</li>
                    </ul>
                </li>
                
                <li><strong>Ringe:</strong> Zahlen werden verwendet, um Verbindungspunkte zu bezeichnen
                    <ul>
                        <li>C1CCCC1 - Cyclopentan</li>
                        <li>c1ccccc1 - Benzol (aromatisch)</li>
                        <li>N1CCOCC1 - Morpholin</li>
                    </ul>
                </li>
                
                <li><strong>Verzweigungen:</strong> Werden durch Klammern bezeichnet
                    <ul>
                        <li>CC(=O)O - Essigsäure</li>
                        <li>CC(C)C - Isobutan</li>
                        <li>NC(=O)C1CC1 - Cyclopropancarboxamid</li>
                    </ul>
                </li>
                
                <li><strong>Stereochemie:</strong> Werden durch spezielle Symbole bezeichnet
                    <ul>
                        <li>C/C=C/C - trans-Buten-2</li>
                        <li>C/C=C\\C - cis-Buten-2</li>
                        <li>C[@H](F)(Cl)Br - chirales Zentrum</li>
                    </ul>
                </li>
            </ol>
            
            <h3 style='color: #3498db;'>Detaillierte Umwandlungsbeispiele:</h3>
            
            <table border='1' cellpadding='8' style='border-collapse: collapse; width: 100%; margin: 15px 0;'>
            <tr style='background-color: #f8f9fa;'>
                <th style='padding: 10px;'>Molekül</th>
                <th style='padding: 10px;'>SMILES</th>
                <th style='padding: 10px;'>Erklärung</th>
            </tr>
            <tr>
                <td style='padding: 8px;'><strong>Wasser</strong></td>
                <td style='padding: 8px;'><code>O</code></td>
                <td style='padding: 8px;'>Einfaches Molekül - nur Sauerstoffsymbol</td>
            </tr>
            <tr>
                <td style='padding: 8px;'><strong>Ethanol</strong></td>
                <td style='padding: 8px;'><code>CCO</code></td>
                <td style='padding: 8px;'>Kette aus 2 Kohlenstoffen und Sauerstoff</td>
            </tr>
            <tr>
                <td style='padding: 8px;'><strong>Benzol</strong></td>
                <td style='padding: 8px;'><code>c1ccccc1</code></td>
                <td style='padding: 8px;'>Aromatischer Ring (Kleinbuchstaben)</td>
            </tr>
            <tr>
                <td style='padding: 8px;'><strong>Aspirin</strong></td>
                <td style='padding: 8px;'><code>CC(=O)OC1=CC=CC=C1C(=O)O</code></td>
                <td style='padding: 8px;'>Komplexe Struktur mit Funktionsgruppen</td>
            </tr>
            <tr>
                <td style='padding: 8px;'><strong>Koffein</strong></td>
                <td style='padding: 8px;'><code>CN1C=NC2=C1C(=O)N(C(=O)N2C)C</code></td>
                <td style='padding: 8px;'>Heterocyclische Ringe mit Methylgruppen</td>
            </tr>
            <tr>
                <td style='padding: 8px;'><strong>Glukose</strong></td>
                <td style='padding: 8px;'><code>C(C1C(C(C(C(O1)O)O)O)O)O</code></td>
                <td style='padding: 8px;'>Komplexe Struktur mit vielen OH-Gruppen</td>
            </tr>
            </table>
            
            <h3 style='color: #3498db;'>Komplexes Beispiel für Parsing:</h3>
            <p><strong>Molekül:</strong> Acetylsalicylsäure (Aspirin)</p>
            <p><strong>SMILES:</strong> <code style='background-color: #f8f9fa; padding: 5px;'>CC(=O)OC1=CC=CC=C1C(=O)O</code></p>
            
            <p><strong>Parsing nach Teilen:</strong></p>
            <ul>
                <li><code>CC</code> - Methylgruppe (CH3-)</li>
                <li><code>(=O)</code> - Carbonylgruppe (C=O)</li>
                <li><code>O</code> - Sauerstoff der Esterbindung</li>
                <li><code>c1ccccc1</code> - Benzolring (aromatisch)</li>
                <li><code>C(=O)O</code> - Carboxylgruppe (-COOH)</li>
            </ul>
            
            <h3 style='color: #3498db;'>Spezielle Symbole und ihre Bedeutungen:</h3>
            <ul>
                <li><strong>[]</strong> - Atome mit spezifischen Eigenschaften: 
                    <a href="https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html">[Na+]</a>, 
                    <a href="https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system#Atoms">[Fe++]</a>, 
                    <a href="https://rdkit.org/docs/GettingStartedInPython.html">[nH]</a>
                </li>
                <li><strong>@</strong> - Stereochemie: 
                    <a href="https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html">@</a>, 
                    <a href="https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system#Stereochemistry">@@</a>
                </li>
                <li><strong>.</strong> - Trennung von Molekülen: 
                    <a href="https://rdkit.org/docs/Cookbook.html">Na.Cl</a> (Salz)
                </li>
                <li><strong>/</strong> und <strong>\\</strong> - stereochemische Bindungen: 
                    <a href="https://www.cheminfo.org/flavor/malaria/Utilities/SMILES_canonization/index.html">C/C=C/C</a>
                </li>
            </ul>
            
            <h3 style='color: #3498db;'>Nützliche Ressourcen und Tools:</h3>
            <ul>
                <li><a href="https://pubchem.ncbi.nlm.nih.gov/">PubChem</a> - große Datenbank chemischer Verbindungen</li>
                <li><a href="https://www.chemspider.com/">ChemSpider</a> - Suche nach chemischen Strukturen</li>
                <li><a href="https://www.rdkit.org/">RDKit</a> - Bibliothek für chemische Informatik</li>
                <li><a href="https://cactus.nci.nih.gov/chemical/structure">NCI Cactus</a> - Konverter von Namen zu SMILES</li>
                <li><a href="https://www.cheminfo.org/flavor/malaria/Utilities/SMILES_canonization/index.html">SMILES Canonizer</a> - Überprüfung und Kanonisierung</li>
                <li><a href="https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system">Wikipedia SMILES</a> - detaillierte Dokumentation</li>
                <li><a href="https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html">Daylight Theory</a> - offizielle Spezifikation</li>
            </ul>
            
            <h3 style='color: #3498db;'>Wie erhalten Sie SMILES für Ihre Moleküle:</h3>
            <ul>
                <li><strong>Chemische Editoren:</strong> 
                    <a href="https://chemaxon.com/products/marvin">MarvinSketch</a>, 
                    <a href="https://www.perkinelmer.com/category/chemdraw">ChemDraw</a>
                </li>
                <li><strong>Online-Tools:</strong> 
                    <a href="https://www.molview.org/">MolView</a>, 
                    <a href="https://www.cheminfo.org/">ChemDoodle</a>
                </li>
                <li><strong>Datenbanken:</strong> 
                    <a href="https://www.drugbank.ca/">DrugBank</a>, 
                    <a href="https://www.ebi.ac.uk/chembl/">ChEMBL</a>
                </li>
                <li><strong>Mobile Apps:</strong> 
                    <a href="https://play.google.com/store/apps/details?id=com.cheminfomatics.molprime">MolPrime+</a>
                </li>
            </ul>
            
            <div style='background-color: #e8f5e8; padding: 15px; border-radius: 5px; margin: 15px 0;'>
            <strong>💡 Tipp:</strong> Für komplexe Moleküle empfehlen wir die Verwendung grafischer Editoren, 
            die automatisch SMILES-Code aus der gezeichneten Struktur generieren.
            </div>
            
            <p style='color: #7f8c8d; font-style: italic;'>
            Diese Hilfe wurde erstellt, um beim Lernen der SMILES-Notation zu helfen. Für ein tieferes Studium 
            empfehlen wir, auf die offizielle Dokumentation und spezialisierte Ressourcen zurückzugreifen.
            </p>
            </div>
            """
        elif lang == 'uk':
            return """
            <div style='font-size: 12px; line-height: 1.6;'>
            <h3 style='color: #3498db;'>Що таке SMILES?</h3>
            <p><strong>SMILES</strong> - це стандартний текстовий формат для представлення хімічних структур, 
            розроблений для зручного введення та обміну хімічною інформацією.</p>
            
            <h3 style='color: #3498db;'>Основні правила синтаксису:</h3>
            <ol>
                <li><strong>Атоми:</strong> Позначаються хімічними символами (C, O, N, Cl, Na тощо)
                    <ul>
                        <li>Великі літери: C, O, N - органічні атоми</li>
                        <li>Маленькі літери: c, o, n - ароматичні атоми</li>
                        <li>Елементи в дужках: [Na+], [Cl-], [Fe++] - іони та спеціальні атоми</li>
                    </ul>
                </li>
                
                <li><strong>Зв'язки:</strong>
                    <ul>
                        <li>Одинарний зв'язок: не вказується або '-' (C-C або CC)</li>
                        <li>Подвійний зв'язок: '=' (C=O - карбонільна група)</li>
                        <li>Потрійний зв'язок: '#' (C#N - нітрильна група)</li>
                        <li>Ароматичні зв'язки: автоматично визначаються для malih літер</li>
                    </ul>
                </li>
                
                <li><strong>Цикли:</strong> Використовуються числа для позначення точок з'єднання
                    <ul>
                        <li>C1CCCC1 - циклопентан</li>
                        <li>c1ccccc1 - бензол (ароматичний)</li>
                        <li>N1CCOCC1 - морфолін</li>
                    </ul>
                </li>
                
                <li><strong>Розгалуження:</strong> Позначаються дужками
                    <ul>
                        <li>CC(=O)O - оцтова кислота</li>
                        <li>CC(C)C - ізобутан</li>
                        <li>NC(=O)C1CC1 - циклопропанакарбоксамід</li>
                    </ul>
                </li>
                
                <li><strong>Стереохімія:</strong> Позначається спеціальними символами
                    <ul>
                        <li>C/C=C/C - транс-бутен-2</li>
                        <li>C/C=C\\C - цис-бутен-2</li>
                        <li>C[@H](F)(Cl)Br - хіральний центр</li>
                    </ul>
                </li>
            </ol>
            
            <h3 style='color: #3498db;'>Детальні приклади перетворення:</h3>
            
            <table border='1' cellpadding='8' style='border-collapse: collapse; width: 100%; margin: 15px 0;'>
            <tr style='background-color: #f8f9fa;'>
                <th style='padding: 10px;'>Молекула</th>
                <th style='padding: 10px;'>SMILES</th>
                <th style='padding: 10px;'>Пояснення</th>
            </tr>
            <tr>
                <td style='padding: 8px;'><strong>Вода</strong></td>
                <td style='padding: 8px;'><code>O</code></td>
                <td style='padding: 8px;'>Проста молекула - просто символ кисню</td>
            </tr>
            <tr>
                <td style='padding: 8px;'><strong>Етанол</strong></td>
                <td style='padding: 8px;'><code>CCO</code></td>
                <td style='padding: 8px;'>Ланцюг з 2 вуглеців та кисню</td>
            </tr>
            <tr>
                <td style='padding: 8px;'><strong>Бензол</strong></td>
                <td style='padding: 8px;'><code>c1ccccc1</code></td>
                <td style='padding: 8px;'>Ароматичне кільце (малі літери)</td>
            </tr>
            <tr>
                <td style='padding: 8px;'><strong>Аспірин</strong></td>
                <td style='padding: 8px;'><code>CC(=O)OC1=CC=CC=C1C(=O)O</code></td>
                <td style='padding: 8px;'>Складена структура з функціональними групами</td>
            </tr>
            <tr>
                <td style='padding: 8px;'><strong>Кафеїн</strong></td>
                <td style='padding: 8px;'><code>CN1C=NC2=C1C(=O)N(C(=O)N2C)C</code></td>
                <td style='padding: 8px;'>Гетероциклічні кільця з метильними групами</td>
            </tr>
            <tr>
                <td style='padding: 8px;'><strong>Глюкоза</strong></td>
                <td style='padding: 8px;'><code>C(C1C(C(C(C(O1)O)O)O)O)O</code></td>
                <td style='padding: 8px;'>Складена структура з багатьма OH-групами</td>
            </tr>
            </table>
            
            <h3 style='color: #3498db;'>Складний приклад розбору:</h3>
            <p><strong>Молекула:</strong> Ацетилсаліцилова кислота (Аспірин)</p>
            <p><strong>SMILES:</strong> <code style='background-color: #f8f9fa; padding: 5px;'>CC(=O)OC1=CC=CC=C1C(=O)O</code></p>
            
            <p><strong>Розбір по частинах:</strong></p>
            <ul>
                <li><code>CC</code> - метильна група (CH3-)</li>
                <li><code>(=O)</code> - карбонільна група (C=O)</li>
                <li><code>O</code> - кисень естерного зв'язку</li>
                <li><code>c1ccccc1</code> - бензольне кільце (ароматичне)</li>
                <li><code>C(=O)O</code> - карбоксильна група (-COOH)</li>
            </ul>
            
            <h3 style='color: #3498db;'>Спеціальні символи та їх значення:</h3>
            <ul>
                <li><strong>[]</strong> - атоми з специфічними властивостями: 
                    <a href="https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html">[Na+]</a>, 
                    <a href="https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system#Atoms">[Fe++]</a>, 
                    <a href="https://rdkit.org/docs/GettingStartedInPython.html">[nH]</a>
                </li>
                <li><strong>@</strong> - стереохімія: 
                    <a href="https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html">@</a>, 
                    <a href="https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system#Stereochemistry">@@</a>
                </li>
                <li><strong>.</strong> - розділення молекул: 
                    <a href="https://rdkit.org/docs/Cookbook.html">Na.Cl</a> (сіль)
                </li>
                <li><strong>/</strong> и <strong>\\</strong> - стереохімічні зв'язки: 
                    <a href="https://www.cheminfo.org/flavor/malaria/Utilities/SMILES_canonization/index.html">C/C=C/C</a>
                </li>
            </ul>
            
            <h3 style='color: #3498db;'>Корисні ресурси та інструменти:</h3>
            <ul>
                <li><a href="https://pubchem.ncbi.nlm.nih.gov/">PubChem</a> - велика база хімічних сполук</li>
                <li><a href="https://www.chemspider.com/">ChemSpider</a> - пошук за хімічними структурами</li>
                <li><a href="https://www.rdkit.org/">RDKit</a> - бібліотека для хімічної інформатики</li>
                <li><a href="https://cactus.nci.nih.gov/chemical/structure">NCI Cactus</a> - конвертер назв у SMILES</li>
                <li><a href="https://www.cheminfo.org/flavor/malaria/Utilities/SMILES_canonization/index.html">SMILES Canonizer</a> - перевірка та канонізація</li>
                <li><a href="https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system">Wikipedia SMILES</a> - детальна документація</li>
                <li><a href="https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html">Daylight Theory</a> - офіційна специфікація</li>
            </ul>
            
            <h3 style='color: #3498db;'>Як отримати SMILES для ваших молекул:</h3>
            <ul>
                <li><strong>Хімічні редактори:</strong> 
                    <a href="https://chemaxon.com/products/marvin">MarvinSketch</a>, 
                    <a href="https://www.perkinelmer.com/category/chemdraw">ChemDraw</a>
                </li>
                <li><strong>Онлайн-інструменти:</strong> 
                    <a href="https://www.molview.org/">MolView</a>, 
                    <a href="https://www.cheminfo.org/">ChemDoodle</a>
                </li>
                <li><strong>Бази даних:</strong> 
                    <a href="https://www.drugbank.ca/">DrugBank</a>, 
                    <a href="https://www.ebi.ac.uk/chembl/">ChEMBL</a>
                </li>
                <li><strong>Мобільні додатки:</strong> 
                    <a href="https://play.google.com/store/apps/details?id=com.cheminfomatics.molprime">MolPrime+</a>
                </li>
            </ul>
            
            <div style='background-color: #e8f5e8; padding: 15px; border-radius: 5px; margin: 15px 0;'>
            <strong>💡 Порада:</strong> Для складних молекул рекомендуємо використовувати графічні редактори, 
            які автоматично генерують SMILES код з намальованої структури.
            </div>
            
            <p style='color: #7f8c8d; font-style: italic;'>
            Ця довідка створена для допомоги у вивченні SMILES нотації. Для глибшого вивчення 
            рекомендуємо звертатися до офіційної документації та спеціалізованих ресурсів.
            </p>
            </div>
            """
        elif lang == 'fr':
            return """
            <div style='font-size: 12px; line-height: 1.6;'>
            <h3 style='color: #3498db;'>Qu'est-ce que SMILES ?</h3>
            <p><strong>SMILES</strong> - est un format texte standard pour représenter les structures chimiques, 
            développé pour une saisie pratique et l'échange d'informations chimiques.</p>
            
            <h3 style='color: #3498db;'>Règles de syntaxe de base :</h3>
            <ol>
                <li><strong>Atomes :</strong> Désignés par des symboles chimiques (C, O, N, Cl, Na etc.)
                    <ul>
                        <li>Lettres majuscules : C, O, N - atomes organiques</li>
                        <li>Lettres minuscules : c, o, n - atomes aromatiques</li>
                        <li>Éléments entre crochets : [Na+], [Cl-], [Fe++] - ions et atomes spéciaux</li>
                    </ul>
                </li>
                
                <li><strong>Liaisons :</strong>
                    <ul>
                        <li>Liaison simple : non indiquée ou '-' (C-C ou CC)</li>
                        <li>Liaison double : '=' (C=O - groupe carbonyle)</li>
                        <li>Liaison triple : '#' (C#N - groupe nitrile)</li>
                        <li>Liaisons aromatiques : automatiquement déterminées pour les lettres minuscules</li>
                    </ul>
                </li>
                
                <li><strong>Cycles :</strong> Des chiffres sont utilisés pour désigner les points de connexion
                    <ul>
                        <li>C1CCCC1 - cyclopentane</li>
                        <li>c1ccccc1 - benzène (aromatique)</li>
                        <li>N1CCOCC1 - morpholine</li>
                    </ul>
                </li>
                
                <li><strong>Branchages :</strong> Désignés par des parenthèses
                    <ul>
                        <li>CC(=O)O - acide acétique</li>
                        <li>CC(C)C - isobutane</li>
                        <li>NC(=O)C1CC1 - cyclopropanecarboxamide</li>
                    </ul>
                </li>
                
                <li><strong>Stéréochimie :</strong> Désignée par des symboles spéciaux
                    <ul>
                        <li>C/C=C/C - trans-butène-2</li>
                        <li>C/C=C\\C - cis-butène-2</li>
                        <li>C[@H](F)(Cl)Br - centre chiral</li>
                    </ul>
                </li>
            </ol>
            
            <h3 style='color: #3498db;'>Exemples détaillés de conversion :</h3>
            
            <table border='1' cellpadding='8' style='border-collapse: collapse; width: 100%; margin: 15px 0;'>
            <tr style='background-color: #f8f9fa;'>
                <th style='padding: 10px;'>Molécule</th>
                <th style='padding: 10px;'>SMILES</th>
                <th style='padding: 10px;'>Explication</th>
            </tr>
            <tr>
                <td style='padding: 8px;'><strong>Eau</strong></td>
                <td style='padding: 8px;'><code>O</code></td>
                <td style='padding: 8px;'>Molécule simple - juste symbole d'oxygène</td>
            </tr>
            <tr>
                <td style='padding: 8px;'><strong>Éthanol</strong></td>
                <td style='padding: 8px;'><code>CCO</code></td>
                <td style='padding: 8px;'>Chaîne de 2 carbones et oxygène</td>
            </tr>
            <tr>
                <td style='padding: 8px;'><strong>Benzène</strong></td>
                <td style='padding: 8px;'><code>c1ccccc1</code></td>
                <td style='padding: 8px;'>Anneau aromatique (lettres minuscules)</td>
            </tr>
            <tr>
                <td style='padding: 8px;'><strong>Aspirine</strong></td>
                <td style='padding: 8px;'><code>CC(=O)OC1=CC=CC=C1C(=O)O</code></td>
                <td style='padding: 8px;'>Structure complexe avec groupes fonctionnels</td>
            </tr>
            <tr>
                <td style='padding: 8px;'><strong>Caféine</strong></td>
                <td style='padding: 8px;'><code>CN1C=NC2=C1C(=O)N(C(=O)N2C)C</code></td>
                <td style='padding: 8px;'>Anneaux hétérocycliques avec groupes méthyle</td>
            </tr>
            <tr>
                <td style='padding: 8px;'><strong>Glucose</strong></td>
                <td style='padding: 8px;'><code>C(C1C(C(C(C(O1)O)O)O)O)O</code></td>
                <td style='padding: 8px;'>Structure complexe avec de nombreux groupes OH</td>
            </tr>
            </table>
            
            <h3 style='color: #3498db;'>Exemple complexe d'analyse :</h3>
            <p><strong>Molécule :</strong> Acide acétylsalicylique (Aspirine)</p>
            <p><strong>SMILES :</strong> <code style='background-color: #f8f9fa; padding: 5px;'>CC(=O)OC1=CC=CC=C1C(=O)O</code></p>
            
            <p><strong>Analyse par parties :</strong></p>
            <ul>
                <li><code>CC</code> - groupe méthyle (CH3-)</li>
                <li><code>(=O)</code> - groupe carbonyle (C=O)</li>
                <li><code>O</code> - oxygène de la liaison ester</li>
                <li><code>c1ccccc1</code> - anneau benzène (aromatique)</li>
                <li><code>C(=O)O</code> - groupe carboxyle (-COOH)</li>
            </ul>
            
            <h3 style='color: #3498db;'>Symboles spéciaux et leurs significations :</h3>
            <ul>
                <li><strong>[]</strong> - atomes avec propriétés spécifiques : 
                    <a href="https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html">[Na+]</a>, 
                    <a href="https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system#Atoms">[Fe++]</a>, 
                    <a href="https://rdkit.org/docs/GettingStartedInPython.html">[nH]</a>
                </li>
                <li><strong>@</strong> - stéréochimie : 
                    <a href="https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html">@</a>, 
                    <a href="https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system#Stereochemistry">@@</a>
                </li>
                <li><strong>.</strong> - séparation des molécules : 
                    <a href="https://rdkit.org/docs/Cookbook.html">Na.Cl</a> (sel)
                </li>
                <li><strong>/</strong> et <strong>\\</strong> - liaisons stéréochimiques : 
                    <a href="https://www.cheminfo.org/flavor/malaria/Utilities/SMILES_canonization/index.html">C/C=C/C</a>
                </li>
            </ul>
            
            <h3 style='color: #3498db;'>Ressources et outils utiles :</h3>
            <ul>
                <li><a href="https://pubchem.ncbi.nlm.nih.gov/">PubChem</a> - grande base de données de composés chimiques</li>
                <li><a href="https://www.chemspider.com/">ChemSpider</a> - recherche par structures chimiques</li>
                <li><a href="https://www.rdkit.org/">RDKit</a> - bibliothèque pour l'informatique chimique</li>
                <li><a href="https://cactus.nci.nih.gov/chemical/structure">NCI Cactus</a> - convertisseur de noms en SMILES</li>
                <li><a href="https://www.cheminfo.org/flavor/malaria/Utilities/SMILES_canonization/index.html">SMILES Canonizer</a> - validation et canonisation</li>
                <li><a href="https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system">Wikipedia SMILES</a> - documentation détaillée</li>
                <li><a href="https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html">Daylight Theory</a> - spécification officielle</li>
            </ul>
            
            <h3 style='color: #3498db;'>Comment obtenir SMILES pour vos molécules :</h3>
            <ul>
                <li><strong>Éditeurs chimiques :</strong> 
                    <a href="https://chemaxon.com/products/marvin">MarvinSketch</a>, 
                    <a href="https://www.perkinelmer.com/category/chemdraw">ChemDraw</a>
                </li>
                <li><strong>Outils en ligne :</strong> 
                    <a href="https://www.molview.org/">MolView</a>, 
                    <a href="https://www.cheminfo.org/">ChemDoodle</a>
                </li>
                <li><strong>Bases de données :</strong> 
                    <a href="https://www.drugbank.ca/">DrugBank</a>, 
                    <a href="https://www.ebi.ac.uk/chembl/">ChEMBL</a>
                </li>
                <li><strong>Applications mobiles :</strong> 
                    <a href="https://play.google.com/store/apps/details?id=com.cheminfomatics.molprime">MolPrime+</a>
                </li>
            </ul>
            
            <div style='background-color: #e8f5e8; padding: 15px; border-radius: 5px; margin: 15px 0;'>
            <strong>💡 Conseil :</strong> Pour les molécules complexes, nous recommandons d'utiliser des éditeurs graphiques, 
            qui génèrent automatiquement le code SMILES à partir de la structure dessinée.
            </div>
            
            <p style='color: #7f8c8d; font-style: italic;'>
            Cette aide est créée pour assister à l'apprentissage de la notation SMILES. Pour une étude plus approfondie, 
            nous recommandons de se référer à la documentation officielle et aux ressources spécialisées.
            </p>
            </div>
            """
