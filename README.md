
# Cellular Automata Framework
```
Dipartimento di Informatica
Corso di Laurea Magistrale in Informatica
```



## Introduzione

L’obiettivo di questo progetto è la progettazione, implementazione e analisi di un fra-
mework per il task Cellular Automata. L’implementazione del modello dovrà essere rea-
lizzata sia in fastflow che con l’utilizzo diretto dei thread c++. Nei requisiti è richiesto
il supporto di una griglia toroidale, in cui l’ultima riga/colonna è connessa con la prima
riga/colonna. La regola per effettuare l’aggiornamento dello stato delle celle si assume sia
la stessa per tutte e non cambi col tempo, dovrà inoltre essere applica simultaneamente
sull’intera griglia. Viene infine aggiunta la possibilità di scrivere su disco una rappresen-
tazione tramite immagine di ogni iterazione, inoltre l’assunzione fatta per questo progetto
è che il tempo di computazione della regola sia bilanciato su tutte le celle.

## Design

Il problema viene riconosciuto come uno stencil pattern dato che sono presenti delle dipen-
denze strutturali, infatti, quando bisogna calcolare lo stato di una cella è necessario leggere
lo stato dei suoi "vicini", il quale a sua volta potrebbe essere sovrascritto. Assumendo
che la regola sia bilanciata su tutte le celle, possiamo semplificare la suddivisione del task
creando dei range nella matrice, che verranno poi assegnati ai workers.



