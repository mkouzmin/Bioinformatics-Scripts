from io import StringIO
from chimerax.ui import HtmlToolInstance
from chimerax.core.tools import ToolInstance
from chimerax.ui.gui import MainToolWindow
from chimerax.core.errors import UserError
from chimerax.ui import MainToolWindow
from Qt.QtCore import Qt, QThread, Signal, Slot
from Qt.QtWidgets import (
    QPushButton, QSizePolicy
    , QVBoxLayout, QHBoxLayout, QComboBox
    , QWidget, QSpinBox, QAbstractSpinBox
    , QStackedWidget, QPlainTextEdit
    , QLineEdit, QLabel, QAbstractItemView
)
from Qt.QtGui import QAction
from .widget_defs import (
    ViewDockResultsTable
    , ViewDockResultsSettings
)

_settings = None

class _BaseTool(ToolInstance):

    SESSION_ENDURING = False
    SESSION_SAVE = True

    help = "help:user/tools/viewdockx.html"

    def __init__(self, session, tool_name,**kw):
        super().__init__(session, tool_name)
        self.structures = []
        #self._html_state = None
        #self._loaded_page = False
        #self._block_updates = False
        self.display_name = "Viewdock - new Version"
        self.tool_window = MainToolWindow(self)
        self.tool_window.fill_context_menu = self.fill_context_menu

    #called after tool is initialized. Sets up data structures.
    def setup(self, structures=None, html_state=None):

        self.session.logger.status("setup called")
        #self._html_state = html_state
        try:
            self._setup(structures)
        except ValueError as e:
            self.delete()
            raise

    #
    # Get list of structures that we are displaying
    #
    def _setup(self, structures):
        self.category_name = None
        self.category_rating = "viewdockx_rating"
        self.category_list = []
        self.structures = []
        session = self.session
        if structures is None:
            # Include structures only if they have viewdock data
            from chimerax.atomic import AtomicStructure
            structures = [s for s in session.models.list(type=AtomicStructure)
                          if hasattr(s, "viewdockx_data") and s.viewdockx_data]
        else:
            structures = [s for s in structures
                          if hasattr(s, "viewdockx_data") and s.viewdockx_data]

        if not structures:
            raise UserError("No suitable models found for ViewDockX")
        self.structures = structures
        # t = session.triggers
        # from chimerax.core.models import REMOVE_MODELS, MODEL_DISPLAY_CHANGED
        # self._remove_handler = t.add_handler(REMOVE_MODELS, self._update_models)
        # self._display_handler = t.add_handler(MODEL_DISPLAY_CHANGED,
        #                                       self._update_display)

        #
        # Make sure every structure has a rating
        #
        for s in structures:
            if self.category_rating not in s.viewdockx_data:
                s.viewdockx_data[self.category_rating] = "3"
            else:
                try:
                    r = int(s.viewdockx_data[self.category_rating])
                    if r < 0 or r > 5:
                        raise ValueError("out of range")
                except ValueError:
                    s.viewdockx_data[self.category_rating] = "3"

        #
        # Get union of categories found in all viewdockx_data attributes
        #
        category_set = set()
        for s in self.structures:
            try:
                category_set.update([key for key in s.viewdockx_data])
            except AttributeError:
                pass
        # "name" and "rating" categories are special cases that we separate out
        name_found = False
        for category in category_set:
            if category.lower() == "name":
                self.category_name = category
                #category_set.remove(category)
                name_found = True
                break

        if name_found != True:
            for s in structures:
                s.viewdockx_data["Name"] = "unnamed"
            category_set.update(["Name"])
            self.category_name = "Name"


        self.category_list = sorted(list(category_set), key=str.lower)
        self.session.logger.status(''.join(self.category_list))


    def delete(self):
        #t = self.session.triggers
        #if self._remove_handler:
        #    t.remove_handler(self._remove_handler)
        #    self._remove_handler = None
        #if self._display_handler:
        #    t.remove_handler(self._display_handler)
        #    self._display_handler = None
        super().delete()

    def _update_models(self, trigger=None, trigger_data=None):
        """ Called to update page with current list of models"""
        if trigger_data is not None:
            self.structures = [s for s in self.structures
                               if s not in trigger_data]
        if not self.structures:
            self.delete()
            return
        #import json
        #columns = self._make_columns()

        #add id column
        self.table.add_column("id", data_fetch = lambda x : str(x.id_string))

        for category in self.category_list:
            #title = self._format_column_title(string)
            self.table.add_column(str(category), data_fetch=lambda x, i=category: x.viewdockx_data.get(i, None))
        #cat_list = self.category_list.insert(0,"id")
        #self.table.sortByColumn(cat_list, Qt.AscendingOrder)
        #js = "%s.update_columns(%s);" % (self.CUSTOM_SCHEME, columns)
        #self.html_view.runJavaScript(js)

        self.table.data = self.structures #updates rows in the qt table.


    def _make_columns(self):
        # Construct separate dictionaries for numeric and text data
        numeric_data = {}
        text_data = {}
        # First make the id and name columns
        id_list = []
        name_list = []
        name_attr = self.category_name
        for s in self.structures:
            id_list.append(s.id_string)
            if name_attr:
                name_list.append(s.viewdockx_data.get(name_attr, "unnamed"))
        text_data["id"] = id_list
        if name_attr:
            text_data["name"] = name_list
        # Now make numeric and text versions for each category
        # If there are more numbers than text, then assume numeric
        for category in self.category_list:
            numeric_list = []
            text_list = []
            num_numeric = 0
            num_text = 0
            for s in self.structures:
                datum = s.viewdockx_data.get(category, None)
                if datum is None:
                    numeric_list.append(None)
                else:
                    try:
                        if int(datum) == datum:
                            numeric_list.append(int(datum))
                        else:
                            numeric_list.append(float(datum))
                        num_numeric += 1
                    except ValueError:
                        numeric_list.append(None)
                        num_text += 1
                text_list.append(datum)
            if num_numeric > num_text:
                numeric_data[category] = numeric_list
            else:
                text_data[category] = text_list
        return { "numeric": numeric_data, "text": text_data }



    def show_toggle(self, model_id):
        on = []
        off = []
        structures = self.get_structures(model_id)
        for s in structures:
            if s in self.structures:
                if s.display:
                    off.append(s)
                else:
                    on.append(s)
        self._show_hide(on, off)

    def show_set(self, model_id, onoff):
        structures = self.get_structures(model_id)
        on = []
        off = []
        for s in structures:
            if s.display != onoff and s in self.structures:
                if onoff:
                    on.append(s)
                else:
                    off.append(s)
        self._show_hide(on, off)


    def show_only(self, model_id):
        on = []
        off = []
        structures = self.get_structures(model_id)
        for s in self.structures:
            onoff = s in structures
            if s.display != onoff:
                if onoff:
                    on.append(s)
                else:
                    off.append(s)
        self._show_hide(on, off)

    def _show_hide(self, on, off):
        if on or off:
            from chimerax.core.commands import concise_model_spec, run
            self._block_updates = True
            cmd = []
            if off:
                models = concise_model_spec(self.session, off)
                cmd.append("hide %s models" % models)
            if on:
                models = concise_model_spec(self.session, on)
                cmd.append("show %s models" % models)
            run(self.session, " ; ".join(cmd))
            self._block_updates = False
            #updates html display self._update_display()





class TableTool(_BaseTool):

    CUSTOM_SCHEME = "vdxtable"
    _name_map = {}

    def __init__(self, session, tool_name, name=None,
                 structures=None, html_state=None):
        if name is None:
            start = 1
            while str(start) in self._name_map:
                start += 1
            name = str(start)
        elif name in self._name_map:
            raise KeyError("ViewDockX name %r already in use" % name)
        self.name = name
        self._name_map[name] = self



        super().__init__(session,"ViewDockX Table (name: %s)" % name)
        #self.setup_page("viewdockx_table.html")

    def setup(self, structures=None, html_state=None):


        super().setup(structures=None, html_state=None)
        #raise UserError("setup called")
        self._build_ui()

    def _make_settings_dict(self, default_cols):
        defaults = {
            title: True for title in default_cols #format titles first
        }
        return defaults

    def _build_ui(self):
        self.tool_window = MainToolWindow(self)
        parent = self.tool_window.ui_area
        global _settings
        if _settings is None:
            _settings = ViewDockResultsSettings(self.session, "Viewdock")
        self.main_layout = QVBoxLayout()
        self.control_widget = QWidget(parent)
        self.buttons_label = QLabel("For chosen entries:", parent=parent)
        self.buttons_widget = QWidget(parent)
        self.button_container = QHBoxLayout()
        self.button_container.addWidget(self.buttons_label)
        self.button_container.addStretch()

        #param_str = self._format_param_str()
        #self.param_report = QLabel(" ".join(["Query:", param_str]), parent)
        self.control_widget.setVisible(False)

        default_col_list = self.category_list
        default_cols = self._make_settings_dict(default_col_list)
        self.table = ViewDockResultsTable(self.control_widget, default_cols, _settings, parent)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self._update_models()
        self.table.launch(suppress_resize=True)
        #if(self._display_selected(self.table.selected)) is None:
        #    raise UserError("display selected result is None")
        self.table.get_selection.connect(lambda: self._display_selected(self.table.selected))

        self.table.selection_changed.connect(lambda: self._display_selected(self.table.selected))

        self.tool_window.fill_context_menu = self.fill_context_menu

        self.graph_button = QPushButton("Graph", parent=self.buttons_widget)
        self.push_button = QPushButton("Plot", parent=self.buttons_widget)
        self.hbonds_button = QPushButton("Add HBonds", parent=self.buttons_widget)
        self.clash_button = QPushButton("Add Clashes", parent=self.buttons_widget)

        #        self.load_button = QPushButton("Load Structures", parent=self.load_buttons_widget)
        self.button_container.addWidget(self.graph_button)
        self.button_container.addWidget(self.push_button)
        self.button_container.addWidget(self.hbonds_button)
        self.button_container.addWidget(self.clash_button)
        self.hbonds_button.clicked.connect(lambda: self._cb_hb(self.table.selected))
        self.clash_button.clicked.connect(lambda: self._cb_clash(self.table.selected))
        #        self.load_db_button.clicked.connect(lambda: self.load_sequence(self.table.selected))

        #self.main_layout.addWidget(self.param_report)
        self.main_layout.addWidget(self.table)
        self.buttons_widget.setLayout(self.button_container)
        self.main_layout.addWidget(self.control_widget)

        #        self.show_best_matching_container = QWidget(parent=parent)
        #        self.show_best_matching_layout = QHBoxLayout()
        #        self.only_best_matching = QCheckBox("List only best chain per PDB", parent=parent)
        #        self.only_best_matching.stateChanged.connect(self._on_best_matching_state_changed)
        #        self.show_best_matching_layout.addStretch()
        #        self.show_best_matching_layout.addWidget(self.only_best_matching)
        #        self.show_best_matching_container.setLayout(self.show_best_matching_layout)
        #        self.show_best_matching_container.setVisible(False)
        #       self.main_layout.addWidget(self.show_best_matching_container)
        self.main_layout.addWidget(self.buttons_widget)

        #        self.save_button_container = QWidget(parent=parent)
        #        self.save_button_layout = QHBoxLayout()
        #        self.save_button = QPushButton("Save Results as TSV")
        #        self.save_button_layout.addStretch()
        #        self.save_button_layout.addWidget(self.save_button)
        #        self.save_button_container.setLayout(self.save_button_layout)
        #        self.save_button.clicked.connect(self.save_as_tsv)
        #        self.main_layout.addWidget(self.save_button_container)

        #for layout in [self.main_layout, self.show_best_matching_layout, self.save_button_layout]:
        for layout in [self.main_layout]:
            layout.setContentsMargins(2, 2, 2, 2)
            layout.setSpacing(2)

        self.tool_window.ui_area.setLayout(self.main_layout)
        self.tool_window.manage('side')

        self.session.logger.status("UI Built")


    def _display_selected(self,selections) -> None:
        ids = [hit.id_string for hit in selections]
        #ids.insert(0,0)
        seqs = []
        on = []
        off = []
        for sid in ids:
            for s in self.structures:
                if s.id_string == sid:
                    seqs.append(s)
        for s in self.structures:
            onoff = s in seqs
            if s.display != onoff:
                if onoff:
                    on.append(s)
                else:
                    off.append(s)
        self._show_hide(on, off)

    #def _update_rows(self): #Called to create a list of lists of sequence data. Not efficient as a dict of dicts of column values. Possibly switch to pandas structured array instead.
        #self.table.data = [[s.viewdockx_data[key] for s in self.structures] for key in self.category_list]
        #self.table.data = {key: viewdockx_data[key] for key in }
        #raise UserError("data chosen")



    def delete(self):
        del self._name_map[self.name]
        super().delete()

    @classmethod
    def find(cls, name):
        if name is None:
            keys = cls._name_map.keys()
            if len(cls._name_map) > 1:
                raise KeyError("ViewDockX name must be specified when "
                               "there are multiple instances.")
            elif len(cls._name_map) == 0:
                raise KeyError("No active ViewDockX instance.")
            return list(cls._name_map.values())[0]
        else:
            try:
                return cls._name_map[name]
            except KeyError:
                raise KeyError("No ViewDockX instance named %s" % name)

    def fill_context_menu(self, menu, x, y):
        seq_action = QAction("Load Structures", menu)
        #seq_action.triggered.connect(lambda: self.load(self.table.selected))
        #seq_view_action.triggered.connect(lambda: self._show_mav(self.table.selected))
        #load_from_db_action.triggered.connect(lambda: self.load_sequence(self.table.selected))
        #save_as_tsv_action.triggered.connect(lambda: self.save_as_tsv())
        menu.addAction(seq_action)
        #menu.addAction(seq_view_action)
        #menu.addAction(load_from_db_action)
        #menu.addAction(save_as_tsv_action)

    def _update_ratings(self, trigger=None, trigger_data=None):
        if trigger_data is None:
            trigger_data = self.structures
        ratings = [(s.atomspec[1:], s.viewdockx_data[self.category_rating])
                   for s in trigger_data]
        #import json
        #js = "%s.update_ratings(%s);" % (self.CUSTOM_SCHEME,
        #                                 json.dumps(rat m

    def handle_scheme(self, url):
        # Called when custom link is clicked.
        # "info" is an instance of QWebEngineUrlRequestInfo
        from urllib.parse import parse_qs
        method = getattr(self, "_cb_" + url.path())
        query = parse_qs(url.query())

                method(query)

    def arrow_down(self):
        self._arrow_key(1)

    def arrow_up(self):
        self._arrow_key(-1)

    def _arrow_key(self, offset):#to do
        js = "%s.arrow_key(%d);" % (self.CUSTOM_SCHEME, offset)
        #self.html_view.runJavaScript(js)

    def _cb_show_all(self, query):
        """shows or hides all structures"""
        self.show_set(None, True)

    def _cb_show_only(self, query):
        """shows only selected structure"""
        try:
            models = query["id"][0]
        except KeyError:
            self.show_set(None, False)
        else:
            self.show_only(models)

    def _cb_rating(self, query):
        """update rating for structure"""
        # May need to fire trigger for notification later
        try:
            model_id = query["id"][0]
            rating = int(query["rating"][0])
        except (KeyError, ValueError):
            return
        structures = self.get_structures(model_id)
        any_change = False
        for s in structures:
            v = str(rating)
            if s.viewdockx_data[self.category_rating] != v:
                s.viewdockx_data[self.category_rating] = v
                any_change = True
        if any_change:
            self._update_ratings(trigger_data=structures)

    def _cb_graph(self, query):
        tool = ChartTool(self.session, "ViewDockX Graph")
        tool.setup(self.structures)

    def _cb_plot(self, query):
        tool = PlotTool(self.session, "ViewDockX Plot")
        tool.setup(self.structures)

    def _cb_hb(self, query):
        self._count_pbonds(query, "hbonds", "hydrogen bonds", "HBonds")

    def _count_pbonds(self, query, finder, cat_name, column_name):
        # Create hydrogen bonds between receptor(s) and ligands
        from chimerax.core.commands import concise_model_spec, run
        from chimerax.atomic import AtomicStructure
        mine = concise_model_spec(self.session, self.structures)
        self.session.logger.status(''.join(mine))
        all = self.session.models.list(type=AtomicStructure)
        others = concise_model_spec(self.session,
                                    set(all) - set(self.structures))
        cmd = ("%s %s restrict %s "
               "reveal true intersubmodel true" % (finder, mine, others))
        run(self.session, cmd)
        self._count_pb(cat_name, column_name)

    def _count_pb(self, group_name, key):
        # Count up the hydrogen bonds for each structure
        pbg = self.session.pb_manager.get_group(group_name)
        pa1, pa2 = pbg.pseudobonds.atoms
        for s in self.structures:
            atoms = s.atoms
            ma1 = pa1.mask(atoms)
            ma2 = pa2.mask(atoms)
            s.viewdockx_data[key] = (ma1 ^ ma2).sum()
        # Make sure HBonds is in our list of columns
        if key not in self.category_list:
            self.category_list.append(key)
            self.category_list.sort(key=str.lower)
        self._update_models()

    def _cb_clash(self, query):
        self._count_pbonds(query, "clashes", "clashes", "Clashes")

    def _cb_export(self, query):
        from chimerax.ui.open_save import SaveDialog
        sd = SaveDialog(self.session, data_formats=[self.session.data_formats["mol2"]])
        if not sd.exec():
            return
        path = sd.get_path()
        if path is None:
            return
        prefix = "##########"
        from chimerax.mol2 import write_mol2
        with open(path, "w") as outf:
            for s in self.structures:
                with OutputCache() as sf:
                    write_mol2(self.session, sf, models=[s])
                for item in s.viewdockx_data.items():
                    print(prefix, "%s: %s\n" % item, end='', file=outf)
                print("\n", end='', file=outf)
                print(sf.saved_output, end='', file=outf)
                print("\n\n", end='', file=outf)

    def _cb_prune(self, query):
        if 'stars' not in query:
            raise UserError("Must select a number of stars next to the 'Close' button")
        stars = int(query["stars"][0])
        structures = [s for s in self.structures
                      if int(s.viewdockx_data[self.category_rating]) <= stars]
        if not structures:
            print("No structures closed")
            return
        self.session.models.close(structures)

    def _cb_columns_updated(self, query):
        #self._update_display()
        self._update_ratings()

    def _cb_arrow(self, query):
        from chimerax.core.commands import run
        direction = query["direction"][0]
        cmd = "viewdock %s name %s" % (direction, self.name)
        run(self.session, cmd, log=False)





class ChartTool(_BaseTool):

    CUSTOM_SCHEME = "vdxchart"

    help = "help:user/tools/viewdockx.html#plots"

    def __init__(self, session, tool_name, structures=None, html_state=None):
        super().__init__(session, "ViewDockX Chart")
        #self.setup_page("viewdockx_chart.html")
        self._build_ui()


    def _build_ui(self):
        self.tool_window = MainToolWindow(self)
        parent = self.tool_window.ui_area
        global _settings
        if _settings is None:
            _settings = ViewDockResultsSettings(self.session, "Viewdock")
        self.main_layout = QVBoxLayout()
        self.control_widget = QWidget(parent)
        self.buttons_label = QLabel("For chosen entries:", parent=parent)
        self.buttons_widget = QWidget(parent)
        self.button_container = QHBoxLayout()
        self.button_container.addWidget(self.buttons_label)
        self.button_container.addStretch()

        #param_str = self._format_param_str()
        #self.param_report = QLabel(" ".join(["Query:", param_str]), parent)
        self.control_widget.setVisible(False)

        default_col_list = self.category_list
        default_cols = self._make_settings_dict(default_col_list)
        self.table = ViewDockResultsTable(self.control_widget, default_cols, _settings, parent)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self._update_models()
        self.table.launch(suppress_resize=True)
        #if(self._display_selected(self.table.selected)) is None:
        #    raise UserError("display selected result is None")
        self.table.get_selection.connect(lambda: self._display_selected(self.table.selected))

        self.table.selection_changed.connect(lambda: self._display_selected(self.table.selected))

        self.tool_window.fill_context_menu = self.fill_context_menu

        self.graph_button = QPushButton("Graph", parent=self.buttons_widget)
        self.push_button = QPushButton("Plot", parent=self.buttons_widget)
        self.hbonds_button = QPushButton("Add HBonds", parent=self.buttons_widget)
        self.clash_button = QPushButton("Add Clashes", parent=self.buttons_widget)

        #        self.load_button = QPushButton("Load Structures", parent=self.load_buttons_widget)
        self.button_container.addWidget(self.graph_button)
        self.button_container.addWidget(self.push_button)
        self.button_container.addWidget(self.hbonds_button)
        self.button_container.addWidget(self.clash_button)
        self.hbonds_button.clicked.connect(lambda: self._cb_hb(self.table.selected))
        self.clash_button.clicked.connect(lambda: self._cb_clash(self.table.selected))
        #        self.load_db_button.clicked.connect(lambda: self.load_sequence(self.table.selected))

        #self.main_layout.addWidget(self.param_report)
        self.main_layout.addWidget(self.table)
        self.buttons_widget.setLayout(self.button_container)
        self.main_layout.addWidget(self.control_widget)

        #        self.show_best_matching_container = QWidget(parent=parent)
        #        self.show_best_matching_layout = QHBoxLayout()
        #        self.only_best_matching = QCheckBox("List only best chain per PDB", parent=parent)
        #        self.only_best_matching.stateChanged.connect(self._on_best_matching_state_changed)
        #        self.show_best_matching_layout.addStretch()
        #        self.show_best_matching_layout.addWidget(self.only_best_matching)
        #        self.show_best_matching_container.setLayout(self.show_best_matching_layout)
        #        self.show_best_matching_container.setVisible(False)
        #       self.main_layout.addWidget(self.show_best_matching_container)
        self.main_layout.addWidget(self.buttons_widget)

        #        self.save_button_container = QWidget(parent=parent)
        #        self.save_button_layout = QHBoxLayout()
        #        self.save_button = QPushButton("Save Results as TSV")
        #        self.save_button_layout.addStretch()
        #        self.save_button_layout.addWidget(self.save_button)
        #        self.save_button_container.setLayout(self.save_button_layout)
        #        self.save_button.clicked.connect(self.save_as_tsv)
        #        self.main_layout.addWidget(self.save_button_container)

        #for layout in [self.main_layout, self.show_best_matching_layout, self.save_button_layout]:
        for layout in [self.main_layout]:
            layout.setContentsMargins(2, 2,  2, 2)
            layout.setSpacing(2)

        self.tool_window.ui_area.setLayout(self.main_layout)

        self.tool_window.manage('side')

        self.session.logger.status("UI Built")

    def handle_scheme(self, url):
        # Called when custom link is clicked.
        # "info" is an instance of QWebEngineUrlRequestInfo
        from urllib.parse import parse_qs
        method = getattr(self, "_cb_" + url.path())
        query = parse_qs(url.query())
        method(query)

    def _cb_show_only(self, query):
        """shows or hides all structures"""
        self.show_only(query["id"][0])

    def _cb_show_toggle(self, query):
        """shows or hides all structures"""
        self.show_toggle(query["id"][0])


class PlotTool(_BaseTool):

    CUSTOM_SCHEME = "vdxplot"

    help = "help:user/tools/viewdockx.html#plots"

    def __init__(self, session, tool_name, structures=None, html_state=None):
        super().__init__(session, "ViewDockX Plot")
        #self.setup_page("viewdockx_plot.html")

    def handle_scheme(self, url):
        # Called when custom link is clicked.
        # "info" is an instance of QWebEngineUrlRequestInfo
        from urllib.parse import parse_qs
        method = getattr(self, "_cb_" + url.path())
        query = parse_qs(url.query())
        method(query)

    def _cb_show_only(self, query):
        """shows or hides all structures"""
        self.show_only(query["id"][0])

    def _cb_show_toggle(self, query):
        """shows or hides all structures"""
        self.show_toggle(query["id"][0])