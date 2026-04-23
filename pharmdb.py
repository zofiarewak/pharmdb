# Duży projekt zaliczeniowy z ASD 2025
# PharmDB – system zarządzania i analizy leków
from collections import deque 
import heapq
class PharmDB:
    # Klasa przechowująca kompleksową bazę danych leków oraz umożliwiająca analizę ich wzajemnych zależności.
    # System zapewnia efektywne odpowiedzi na zapytania dotyczące bezpieczeństwa przyjmowanych leków 
    # i optymalizacji ich doboru.
    #
    # Kazdy lek jest opisany następującymi atrybutami:
    # - id: unikalny identyfikator nadawany przez system (np. "D0001", "D0002")
    # - name: unikalna nazwa leku (np. "Apap", "Ibuprom")
    # - indications: lista wskazań terapeutycznych (chorób) wraz z poziomem skuteczności (skala 1-10)
    # - substitutes: lista leków, które mogą być zastąpione przez dany lek
    #               (jeśli lek A ma na liście substitutes lek B, oznacza to, że B może być zastąpiony przez A)
    # - side_effects: lista działań niepożądanych (nazwa objawu, poziom dolegliwości 1-3, częstotliwość)

    def __init__(self):
        # Inicjalizuje bazę danych leków oraz odpowiednie struktury danych.
        # Wymagana złożoność czasowa: O(1)

        self.baza = {} #tu przechowuje wszystkie id leków, aby wiedzieć jakie leki mamy w naszej bazie
        self.indications_heap = {} #slownik, ktory przechowuje nam heap dla kazdej choroby danego leku. Przechowujemy tam informację o skuteczności i indeksie drug_id, oraz time_added (kiedy lek został dodany do bazy). 
        self.efficacy_number = {} #slownik ktory dla danego drug_id przechowuje tablice z informacja ile jest wskazan dla danego efficacy   
        self.alternative_drugs = {} #dla kazdego leku przechowuje zbiór lekow ktore moga go zastapic 
        self.side_effects_dict = {} #dla każdego leku przechowuje najgorszy skutek uboczny, jego risk_score i time_added (by pamiętać który lek został dodany wcześniej)
        
        self.avl_root = None #Przechowuję korzeń, który trzyma niezbędne informacje do dalszego korzystania z AVL trees. Ma informację o lewym i prawym dziecku, więc ostatecznie stworzy się nam drzewo AVL. 

        self.next_number = 1 #inicjalizacja ostatniej liczby w indeksie leku- pierwszy to A0000, bedziemy to zwiekszac w add_drug 
        self.last_letter = 'A' #analogicznie litera-> zaczynamy od A. 
        self.time_added = 0

        #Elementy, które potrzebne mi są do longest_alternative_list, aby uzyskać zadaną złożoność:
        self.longest_path_drug_id = {} #Słownik, który dla każdego leku przechowywać będzie swojego poprzednika w najdłuższej ścieżce danego leku, oraz jej długość (w postaci krotki)
        self.longest_path_length = 0 #Tu przechowuję informację o tym, jaka jest dotychczasowa długość najdłuższej ścieżki alternatyw. 
        self.previous_drug_in_longest_list = None #Tutaj natomiast mam informację o tym, na czym się ta najdłuższa lista kończy. 

    def add_drug(self, drug_name, indications=None, substitutes=None, side_effects=None):

        '''Funkcja odpowiada za dodawanie leków do bazy danych. Wypełnia listę id w bazie, a także słowniki potrzebne w kolejnych funkcjach.
        Zwraca ona id dodanego leku w postaci litery, oraz 4 cyfr. Id przydzielane są w porządku leksykograficznym.
        '''

        # Dodaje nowy lek do bazy danych o podanych wskazaniach terapeutycznych, 
        # lekach zastępowanych i efektach ubocznych.
        #
        # Args:
        #   drug_name (str): nazwa leku
        #   indications (list, optional): lista krotek (choroba, skuteczność)
        #   substitutes (list, optional): lista identyfikatorów leków, które mogą być zastąpione przez nowy lek
        #   side_effects (list, optional): lista krotek (objaw, poziom dolegliwości, częstotliwość)
        #
        # Returns:
        #   str: identyfikator dodanego leku
        #
        # Wymagana złożoność czasowa: O(k log K + s + e), gdzie:
        #   - k to liczba wskazań terapeutycznych
        #   - K to maksymalna liczba leków dla dowolnego wskazania terapeutycznego
        #   - s to liczba zamienników
        #   - e to liczba działań niepożądanych

        #Jeśli któraś z list to None, to zastępuję pustą listą:
        indications = indications or []
        substitutes = substitutes or []
        side_effects = side_effects or []

        #Identyfikator tworzę jako pierwszy człon dając zmienną last_letter, a później next_number, ale czterocyfrową-> puste miejsca zostaną wypełnione zerami. 
        drug_id = f"{self.last_letter}{self.next_number:04d}"
        drug = [drug_id, drug_name, indications,substitutes,side_effects]
        self.baza[drug_id] = drug #zapisuje lek w bazie
        
        if self.next_number >= 10000: #Jesli nasza liczba przekroczy limit, przechodzimy do kolejnej litery
            if self.last_letter !='Z': #Czy mieścimy się a alfabecie? (Czy limit indeksów nie został przekroczony)
                self.last_letter = chr(ord(self.last_letter) + 1) #Zwiekszamy literkę (zamieniamy na liczbe, dodajemy 1 i ponownie zamieniamy na litere)
                self.next_number = 0 
            else: raise Exception("Zbyt dużo indeksów w bazie!")

        self.next_number += 1 #Zwiększamy next_number o 1
        
        efficacy_number_dict = self.efficacy_number
        efficacy_number_dict[drug_id] = [0] * 11 #Inicjalizacja tablicy w słowniku
        
        current_time = self.time_added
        self.time_added += 1

        #Tworzenie heapq (kopca minimalnego) -> zamienię sobie wartosc efficacy na ujemna, aby móc trzymac najlepsze leki na górze.  
        for d, e in indications: 
            heapq.heappush(self.indications_heap.setdefault(d, []),(-e, -(current_time), drug_id)) #Przechowuje informacje w kolejnosci -efficacy, -czas, drug_id, zeby w razie rownej efficacy kopiec porownywal czasy
            
            for i in range(e + 1): #dla kazdego i w tablicy przechowujemy info o tym ile mamy wskazan o min takiej efficacy, jesli przeszukujemy dana efficacy to musimy zwiekszyc kazdy indeks i (0, efficacy) o 1, bo jesli mamy takie efficacy, to kazdy indeks po drodze zawiera sie w danym efficacy. 
                efficacy_number_dict[drug_id][i] += 1 

        max_path_length = 0 #Jaka jest długość najdłuższej ścieżki kończącej się na tym leku
        previous_drug = None #Jaki jest poprzedzający lek?

        for s in substitutes: #Przeszukujemy leki, które nasz lek może zastąpić

            if s not in self.baza: raise Exception("Zastępowany lek musi być w bazie!")
            
            self.alternative_drugs.setdefault(s, set()).add(drug_id) #Zapisuje drug_id w secie w słowniku alternative_drugs korzystając z setdefault()

            if s in self.longest_path_drug_id: 
                length, _ = self.longest_path_drug_id[s]
                if length > max_path_length:
                    max_path_length = length
                    previous_drug = s

        #Aby zachować jak najlepszą złożoność muszę pamiętać ostatni element obecnej najdłuższej ścieżki i jej długość
        if max_path_length + 1 > self.longest_path_length: 
            self.previous_drug_in_longest_list = drug_id #Zapisuję sobie ostatni element najdłuższej ścieżki 
            self.longest_path_length = max_path_length + 1 #Zapisuję również długość obecnej najdłuższej ścieżki
        
        self.longest_path_drug_id[drug_id] = (max_path_length + 1, previous_drug) #Dodajemy do słownika nasz drug_id, zwiększamy naszą max długość o 1 i zapisujemy previous drug. 

        #Rozważam skutki uboczne:
        if drug_id not in self.side_effects_dict:
            self.side_effects_dict[drug_id] = (None, 0, 0) #Inicjalizacja krotki 
        
        worst, risk_score, _ = self.side_effects_dict[drug_id]

        #Wartości najwyższego poziomu i częstotliwosci skutku:
        worst_p = 0
        worst_c = 0 
        for o, p, c in side_effects:
            
            risk_score += (p*c)
            if (p > worst_p) or (p == worst_p and c > worst_c): #Sprawdzamy czy nasze wartości są wyższe niż obecne najgorsze. Jeśli tak, nadpisujemy zmienne worst_p i worst_c 
                worst = o 
                worst_c = c 
                worst_p = p 

        self.side_effects_dict[drug_id] = (worst, risk_score, current_time)

        #Dla każdego objawu nieporządanego (jeśli nie ma jeszcze takiej częstotliwości) w side_effects dodaję nowy węzeł, korzystając z metody insert_avl. Umożliwia mi to uzyskanie odpowiedniej złożoności
        for side_effect, _, freq in side_effects:
            self.avl_root = self.insert_avl(self.avl_root, freq, drug_id, side_effect)
        return drug_id

    #Funkcje pomocnicze dla mojego słowniko-drzewa:
    
    #Obliczanie wysokości węzła:
    def height(self, node):    
        if node is None:        
            return 0    
        return node[4] 
    #Balance factor węzła:
    def balance_factor(self, node):
        if node is None:
            return 0
        return self.height(node[2]) - self.height(node[3])

    #Definiuję też rotację prawą i lewą- niezbędne w balansowaniu drzewa AVL. Oddzielne funkcje umożliwią mi działanie rekurencyjne w późniejszym kodzie. 
    def left_rotation(self, node):
        right_node = node[3] 
        node[3] = right_node[2] #right_node[2] to lewe dziecko prawego korzenia, prawe dziecko nim się staje
        right_node[2] = node #obecny węzeł to lewe dziecko swojego dawnego prawego dziecka 

        #Muszę też zaktualizować wysokość:
        node[4] = 1 + max(self.height(node[2]), self.height(node[3]))
        right_node[4] = 1 + max(self.height(right_node[2]), self.height(right_node[3]))

        return right_node 

    def right_rotation(self, node):
        #Funkcja działa analogicznie do left_rotation
        left_node = node[2]
        node[2] = left_node[3]
        left_node[3] = node

        node[4] = 1 + max(self.height(node[2]), self.height(node[3]))
        left_node[4] = 1 + max(self.height(left_node[2]), self.height(left_node[3]))

        return left_node


        
    def balancing_avl(self, node):
        '''
        Funkcja odpowiada za balansowanie drzewa. Pobieram balance factor dla rozpatrywanego węzła. 
        Jeśli nasz bf jest > 1 musze przeprowadzić prawą rotację, natomiast jeśli jest <-1, muszę przeprowadzić lewą rotację.
        Rozważam też sytuację, w której dziecko jest przechylone w przeciwną stronę niż rodzic. Wtedy muszę przeprowadzić podwójną rotację. 
        To działanie wywodzi się z definicji drzew AVL, gdzie balance factor musi należeć do przedziału [-1, 1]. 
        '''
        bf = self.balance_factor(node)

        if bf > 1:
            #Muszę też sprawdzić, czy wysokość prawego poddrzewa nie jest większa niż lewego (warunek < 0 z definicji balance factor).
            #Jeśli tak, to przeprowadzam lewą rotację lewego dziecka aktualnego węzła. 
            if self.balance_factor(node[2]) < 0:
                node[2] = self.left_rotation(node[2])
            return self.right_rotation(node) #Po takim działaniu mogę bezpiecznie przeprowadzić balansowanie drzewa.

        if bf < -1:
            #Tutaj działam analogicznie do sytuacji powyżej.
            if self.balance_factor(node[3]) > 0:
                node[3] = self.right_rotation(node[3])
            return self.left_rotation(node)

        return node #Zwracam zbalansowany węzeł


    def insert_avl(self, node, freq, drug_id, side_effect):
        '''
        Funkcja odpowiada za dodawanie węzła do drzewa avl. 
        '''
        if node is None:
            #Tu tworzę nowy węzeł. Zwracam strukturę węzła jaką ustaliłam, czyli [częstotliwość, [lista nazw leków i objawów dla takiej częstotliwości, sąsiadów (tu nowy węzeł-> brak sąsiadów), wysokość]
            return [freq, [(self.baza[drug_id][1], side_effect)], None, None, 1]

        #W tej części kodu zapisuję lewe i prawe dziecko węzła, w zależności od tego w którą stronę musimy iść (czy dodawane freq jest mniejsze czy większe niż poprzednik)
        if freq < node[0]:
            node[2] = self.insert_avl(node[2], freq, drug_id, side_effect)
        elif freq > node[0]:
            node[3] = self.insert_avl(node[3], freq, drug_id, side_effect)
        else:
            #Jeśli taka częstotliwość jest już w drzewie, to dodaję nazwę leku i objaw nieporządany do listy krotek dla takiej częstotliwości. 
            node[1].append((self.baza[drug_id][1], side_effect))
            return node

        #Obliczam wysokość na podstawie dzieci węzła:
        node[4] = 1 + max(self.height(node[2]), self.height(node[3]))
        return self.balancing_avl(node) #Zwracam zbalansowany węzeł


    def number_of_indications(self, drug_id, min_efficacy):

        '''Funkcja korzysta ze słownika efficacy_number, który przechowuje liczbę wskazań terapeutycznych. 
        Dla każdego leku słownik ten przechowuje tablicę, która dla i-tego pola zawiera informację o liczbie wskazań o efficacy i.
        Takie działanie umożliwia złożoność O(1), gdyż w funkcji pobieramy tylko informację ze słownika. 
        Zwraca ona zmienną indications_number, czyli wartość ze słownika dla danego drug_id, i danego min_efficacy.
        '''

        # Zwraca liczbę wskazań terapeutycznych o efektywności co najmniej min_efficacy dla podanego leku.
        #
        # Args:
        #   drug_id (str): identyfikator leku
        #   min_efficacy (int): minimalna wymagana efektywność
        #
        # Returns:
        #   int: liczba wskazań terapeutycznych spełniających kryterium
        #
        # Wymagana złożoność czasowa: O(1)
        
        #Dzięki słownikowi efficacy_number możemy wywołać liczbę takich wskazań ze złożonością O(1):
        indications_number = self.efficacy_number[drug_id][min_efficacy] 

        return indications_number

    def number_of_alternative_drugs(self, drug_id):
        # Zwraca liczbę leków, które mogą bezpośrednio zastąpić dany lek.
        #
        # Args:
        #   drug_id (str): identyfikator leku
        #
        # Returns:
        #   int: liczba leków mogących zastąpić dany lek
        #
        # Wymagana złożoność czasowa: O(1)

        '''
        Funkcja korzysta ze słownika alternative_drugs, który dla każdego leku przechowuje zbiór jego zastępników. 
        Sprawdzana jest długość zbioru i zwracana w postaci zmiennej alternative_drugs_number. Jeśli dany lek nie występuje w słowniku zwracane jest 0. 
        '''

        if drug_id in self.alternative_drugs:
            alternative_drugs_number = len(self.alternative_drugs[drug_id])
            return alternative_drugs_number
        else: return 0 # w przeciwnym wypadku nie ma jak zastąpić leku-> zwracamy 0. 
        
    def worst_side_effect(self, drug_id):
        '''
        Korzystając ze słownika side_effects_dict, który dla każdego leku przechowuje krotkę, która poza risk_score i time_added zawieraja najgorszy efekt. 
        Zwracamy zmienną worst_effect jeżeli słownik nie jest pusty i podane drug_id jest w słowniku. W przeciwnym wypadku zwracane jest None
        '''

        # Zwraca nazwę najbardziej dotkliwego skutku ubocznego dla podanego leku.
        # Najbardziej dotkliwy to skutek o największej częstotliwości spośród tych o największym poziomie dolegliwości.
        #
        # Args:
        #   drug_id (str): identyfikator leku
        #
        # Returns:
        #   str: nazwa najbardziej dotkliwego skutku ubocznego lub None jeśli lek nie ma skutków ubocznych
        #
        # Wymagana złożoność czasowa: O(1)
        #Zwracamy ten efekt uboczny, który ma najwyższy poziom dolegliwości
        if drug_id in self.side_effects_dict and self.side_effects_dict[drug_id]:
            worst_effect, _, _ = self.side_effects_dict[drug_id]
            return worst_effect
        else: return None

    def risk_score(self, drug_id):
        # Zwraca wskaźnik ryzyka (risk score) zdefiniowany jako ważona suma wszystkich działań niepożądanych leku.
        # Risk score = suma(częstość występowania × poziom dolegliwości) po wszystkich efektach ubocznych.
        #
        # Args:
        #   drug_id (str): identyfikator leku
        #
        # Returns:
        #   float: wskaźnik ryzyka
        #
        # Wymagana złożoność czasowa: O(1)

        '''
        Korzystając z tego samego słownika, co w funkcji worst_side_effect(), tym razem z krotki pobieramy risk_score, czyli policzony w funkcji add_drug wskaźnik ryzyka.
        Zwracamy go, dzięki słownikowi tu również została uzyskana złożoność stała, gdyż jedyne co robi ta funkcja, to pobranie danych z side_effects_dict. 
        '''
        if drug_id in self.side_effects_dict:
            _, risk_score, _= self.side_effects_dict[drug_id]
            
            return risk_score

        else: return 0


    def find_best_alternative(self, drug_id, max_steps=2):
        # Zwraca identyfikator leku o minimalnym ryzyku spośród leków, które można zastosować 
        # zamiast leku drug_id przy ograniczeniu liczby zamian (max_steps).
        #
        # Args:
        #   drug_id (str): identyfikator leku
        #   max_steps (int, optional): maksymalna liczba zamian, domyślnie 2
        #
        # Returns:
        #   str: identyfikator leku o minimalnym ryzyku
        #
        # Wymagana złożoność czasowa: funkcja powinna działać istotnie szybciej niż O(D+S)
        
        '''
        W tej funkcji korzystając z algorytmu BFS, przeszukujemy listę alternatywnych leków dla drug_id pobraną z funkcji alternative_drugs.
        Listę leków traktujemy jak graf, który przeszukujemy wszerz.
        Substitutes to sąsiedzi drug_id, który jest source node w algorytmie. 
        - visited to zbiór, który zapamiętuje wszystkie odwiedzone już leki. 
        - todo to deque zawierająca krotki (drug_id, steps). Zawiera obecnie najlepszy lek, oraz liczbę kroków, której również musimy pilnować, by nie przekroczyć max_steps. 
        - min_risk inicjalizujemy jako nieskończoność i z każdym kolejnym przeszukaniem, aktualizujemy wartość. 
        - best_drug to zmienna, która przechowuje obecny najlepszy lek, pod względem risk_score, pobieranym z side_effects_dict. 
        - first_added - przechowuje informację o obecnym najwcześniej dodanym leku. 
        Pobieramy też time_added, aby w razie tego samego risk_score, zwrócić wcześniej dodany lek
        Funkcja zwraca zmienną best_drug. 
        '''

        visited = set()
        todo = deque() 
        min_risk = float("inf")
        best_drug = drug_id #Inicjalizuję nasz best_drug, będzie przechowywał obecny najlepszy lek
        visited.add(drug_id) #ten juz odwiedzilismy, to nasz "source node"

        todo.append((best_drug, 0)) #Kroki narazie to 0
        first_added = float('inf') #Inicjalizacja zmiennej first_added
        
        while todo:

            drug, steps = todo.popleft() 

            if drug in self.side_effects_dict:
                _, risk, time_added = self.side_effects_dict[drug]
                #Opcja 1- kolejny lek różni się ryzykiem
                if risk < min_risk: 
                    min_risk = risk 
                    best_drug = drug
                    first_added = time_added

                #Opcja 2- to samo ryzyko, sprawdzamy czas dodania i wybieramy lek dodany wcześniej jako best_drug
                elif risk == min_risk and time_added < first_added:
                    first_added = time_added
                    best_drug = drug

                #Sprawdzamy czy mieścimy się w krokach do wykonania, po czym iterujemy substitutes dla best_drug, czy te węzły zostały już odwiedzone, jeśli nie to dodaję je do visited. 
                if steps < max_steps: 
                    if drug in self.alternative_drugs:
                        for substitute in self.alternative_drugs[drug]: 
                            if substitute not in visited:
                                visited.add(substitute)
                                todo.append((substitute, steps + 1))
        
        return best_drug

    def longest_alternative_list(self):
        # Zwraca listę identyfikatorów leków stanowiącą najdłuższy ciąg zamienników leków,
        # gdzie każdy kolejny lek na liście może bezpośrednio zamienić lek poprzedni.
        #
        # Przykład: jeśli A może być zastąpione przez B, B przez C lub D, a C przez D,
        # to wynikiem działania tej funkcji powinna być lista [id(A), id(B), id(C), id(D)].
        #
        # Returns:
        #   list: lista identyfikatorów leków
        #
        # Wymagana złożoność czasowa: O(d), gdzie d to długość zwracanej listy

        '''
        Funkcja korzysta z:
        - previous_drug_in_longest_list (zmienna przechowująca ostatni lek z obecnej najdłuższej listy)
        - longest_path_drug_id (słownik który dla każdego drug_id trzyma (długość ścieżki, poprzednik leku))
        Pobieramy ostatni lek z najdłuższej ścieżki, odbudowujemy ją od końca. Za pomocą while dla każdego badanego drug dodajemy go do odbudowywanej listy longest_alternative_list_.
        Jeśli drug jest w słowniku, to pobieramy poprzednika, czyli w naszej liście kolejny lek. Przypisujemy wartość drug = previous_drug. 
        Na koniec z użyciem funkcji reverse() odwracamy listę.
        Zwracamy longest_alternative_list_. 
        '''

        longest_alternative_list_ = []
        drug = self.previous_drug_in_longest_list #Rozważamy kolejno końce naszej ścieżki leków-> na koniec będziemy więc musieli odwrócić naszą listę

        while isinstance(drug, str) and drug: #Czy drug jest w postaci stringa? Upewniam się, czy nie jest none.
            longest_alternative_list_.append(drug)
            if drug in self.longest_path_drug_id:
                _, previous_drug = self.longest_path_drug_id[drug]
            else: 
                previous_drug = None

            drug = previous_drug #Rozważamy leki od końca, więc teraz lek jest poprzednim. 

        longest_alternative_list_.reverse()
        return longest_alternative_list_

    def find_best_drug_for_indication(self, disease_name):
        # Zwraca identyfikator leku o największej efektywności dla wskazanej choroby.
        # W przypadku wielu leków o takiej samej efektywności, zwraca najpóźniej dodany do bazy.
        #
        # Args:
        #   disease_name (str): nazwa choroby
        #
        # Returns:
        #   str: identyfikator leku o największej efektywności
        #
        # Wymagana złożoność czasowa: O(1)

        '''
        Korzystając z kopca minimalnego utworzonego w add_drug, jesteśmy w stanie pobrać wartość best_drug z góry kopca. 
        Dzięki zapisywaniu wartości efficacy na minusie, najlepszy lek jest na górze kopca
        '''

        if disease_name not in self.indications_heap or not self.indications_heap[disease_name]: 
            return None
            
        best_drug = self.indications_heap[disease_name][0][2]

        return best_drug

    def update_best_indication(self, disease_name, new_efficacy):
        # Zmienia efektywność najlepszego leku dla wskazanej choroby, tj.
        # leku, który jest wynikiem działania find_best_drug_for_indication
        #
        # Args:
        #   disease_name (str): nazwa choroby
        #   new_efficacy (int): nowa efektywność
        #
        # Wymagana złożoność czasowa: O(log K)

        ''' 
        Funkcja aktualizuje wartość efficacy najlepszego leku dla danej choroby. 
        Korzystając z funkcji heapreplace() zapisujemy nową wartość efektywności na minusie, by na górze kopca minimalnego znajdywały się leki o najwyższej efektywności. 
        Muszę też odświeżyć słownik efficacy_number, dlatego jeśli nowa efficacy jest większa niż poprzednia, dodajemy 1 w tablicy, aż do indeksu odpowiadającego za new_efficacy. 
        W razie nowej efektywności mniejszej niż poprzednia, alternatywnie odejmujemy 1 od każdego pola tablicy. Pomoże to uniknąć błędów w funkcji number_of_indications
        '''

        best_drug = self.indications_heap[disease_name][0][2]
        
        #Aktualizuję tablicę, przez to, że max efficacy to 10, to zakres będzie max 11 (stały, dzięki czemu mamy wymaganą złożoność)
        old_efficacy = -self.indications_heap[disease_name][0][0] 

        if new_efficacy > old_efficacy:
            for i in range(old_efficacy + 1, new_efficacy + 1):
                self.efficacy_number[best_drug][i] += 1
        else:
            for i in range(new_efficacy + 1, old_efficacy + 1):
                self.efficacy_number[best_drug][i] -= 1

        #Update kopca zawierającego informację o indications dla danego leku 
        heapq.heapreplace(self.indications_heap[disease_name], (-new_efficacy, self.time_added, best_drug)) 
    
    def count_drugs_with_side_effect_frequency(self, min_f, max_f):
        # Zwraca liczbę par (lek, objaw nieporządany) w bazie danych, gdzie lek
        # powoduje objaw niporządany we wskazanym (obustronnie domkniętym) zakresie częstotliwości występowania.
        # Funkcja powinna działać w czasie zamortyzowanym O(log F),
        # gdzie F to sumaryczna liczba działań niepożądanych dla wszystkich leków w bazie danych.

        '''
        Korzystam z własności AVL tree i definiując funkcję pomocniczą count_drugs_with_side_effect_frequency_ używam rekurencji i zliczam węzły z odpowiednim freq dla każdego poddrzewa, zaczynając od korzenia. 
        '''
        def count_drugs_with_side_effect_frequency_(node): 
            if not node:
                return 0

            node_freq = node[0]  #Pobieram wartość częstotliwości dla węzła
            
            if node_freq < min_f: #Idę w prawo (większe wartości freq niż rozpatrywany node), lub lewo (mniejsze wartości), w zależności od tego czy node_freq jest mniejsze, czy wieksze. 
                return count_drugs_with_side_effect_frequency_(node[3])  
            elif node_freq > max_f:
                return count_drugs_with_side_effect_frequency_(node[2])
            else: #Jeśli się zgadza, zliczamy lewe i prawe freq
                left = count_drugs_with_side_effect_frequency_(node[2]) 
                right = count_drugs_with_side_effect_frequency_(node[3])  #prawe dziecko
                return 1 + left + right #Obecny węzeł i lewe i prawe 

        return count_drugs_with_side_effect_frequency_(self.avl_root) #Funkcję wywołuję od korzenia drzewa. 

    def list_drugs_with_side_effect_frequency(self, min_freq, max_freq):
        # Zwraca listę par (lek, objaw niepożądany), dla których częstotliwość występowania objawu
        # mieści się we wskazanym zakresie. Funkcja powinna działać w czasie zamortyzowanym O(log F + m),
        # gdzie m jest liczbą par (lek, objaw) z zadanego przedziału.
        '''
        Analogicznie do funkcji powyżej, korzystam z własności AVL tree i definiując funkcję pomocniczą list_drugs_with_side_effect_frequency_ rozszerzam listę par dla pasujących do przedziału częstotliwości. 
        '''
        def list_drugs_with_side_effect_frequency_(node): #Potrzebuję rekurencji, tak jak w przykładzie powyżej muszę sprawdzać dziecko lewe i prawe w drzewie dla każdego badanego węzła. 
            if not node:
                return []

            node_freq = node[0]

            #W zależności od tego, czy freq rozpatrywanego węzła jest mniejsza czy większa, będę iść w stronę lewą (mniejsze wartości) lub prawą (większe wartości) drzewa AVL. 
            if node_freq < min_freq:
                return list_drugs_with_side_effect_frequency_(node[3])  #idź w prawo, bo tam freq większe.
            elif node_freq > max_freq:
                return list_drugs_with_side_effect_frequency_(node[2])  #idź w lewo, bo tam freq mniejsze, analogicznie, taka jest struktura AVL. 
            else: #Jeśli się zgadza, zliczamy lewe i prawe freq
                pairs = [] #Lokalna lista par dla węzła, którą rozszerzam o listy dla lewego i prawego dziecka
                left_pairs = list_drugs_with_side_effect_frequency_(node[2]) 
                right_pairs = list_drugs_with_side_effect_frequency_(node[3]) 
                pairs.extend(node[1]) #Dodaję obecną parę węzła
                pairs.extend(left_pairs) #Muszę rozszerzyć listę o pairs dla prawego i lewego dziecka 
                pairs.extend(right_pairs)
                return pairs
        return list_drugs_with_side_effect_frequency_(self.avl_root) #Wywołuję funkcję dla korzenia drzewa. 


'''
Źródła:
- Zrozumienie drzew AVL- https://eduinf.waw.pl/inf/alg/001_search/0119.php, filmik na youtube: " AVL Trees Simply Explained " z kanału Maaneth De Silva
- Algorytm tych drzew i modyfikacja pod przykład: insercja do drzewa (w tym balansowanie, right, left rotation i balance_factor)- https://www.geeksforgeeks.org/insertion-in-an-avl-tree/
''' 