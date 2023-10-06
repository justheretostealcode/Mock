class Plasmid(list):
    def __init__(self):
        [self.append(ChromosomalLocation(iX)) for iX in range(2)]

        self.chromosomal_locations = self

    def to_json(self):
        json_rep = [elem.to_json() for elem in self]
        return json_rep

    @classmethod
    def from_json(cls, json_data):
        plasmid = cls()
        for i, elem_data in enumerate(json_data):
            plasmid[i] = ChromosomalLocation.from_json(elem_data, id=i)
        return plasmid

    def __str__(self):
        str_rep = ""
        for elem in self:
            str_rep += str(elem)
            str_rep += "\n"
        return str_rep


class ChromosomalLocation(list):
    def __init__(self, id):
        self.id = id
        # [self.append(Position(iX)) for iX in range(16)]
        [self.append(None) for iX in range(8)]
        self.positions = list(self)

        pass

    def to_json(self):
        json_rep = [elem.to_json() if elem is not None else None for elem in self]
        return json_rep

    @classmethod
    def from_json(cls, json_data, id):
        location = cls(id)
        for i, elem_data in enumerate(json_data):
            location[i] = None
            if elem_data is not None:
                location[i] = Position.from_json(elem_data, id=i)

        return location

    def __str__(self):
        str_rep = "|"
        for elem in self:
            str_rep += " "
            str_rep += str(elem)
            str_rep += " |"
        return str_rep


class Position(dict):
    def __init__(self, id, promoter=None, coding_sequence=None):
        self.id = id

        self.promoter = promoter
        self.coding_sequence = coding_sequence

        self.update(self.to_json())
        pass

    def to_json(self):
        json_dict = {"Promoter": self.promoter,
                     "Coding Sequence": self.coding_sequence}
        return json_dict

    @classmethod
    def from_json(cls, json_data, id):
        return cls(
            id=id,
            promoter=json_data["Promoter"],
            coding_sequence=json_data["Coding Sequence"]
        )

    def __str__(self):
        str_rep = ""
        str_rep += "Promoter: "
        str_rep += str(self.promoter)
        str_rep += ", CDS: "
        str_rep += str(self.coding_sequence)

        return str_rep
