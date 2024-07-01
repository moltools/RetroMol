# -*- coding: utf-8 -*-

"""The nodes module."""

import re
import typing as ty
from abc import ABC, abstractmethod

from retromol.npkg.connection import Neo4jConnection


class Node(ABC):
    """The Node abstract base class."""


class BioactivityLabel(Node):
    """The BioactivityLabel node."""

    @classmethod
    def set_constraints(cls, conn: Neo4jConnection) -> None:
        """Set constraints on the Neo4j database.

        :param conn: The Neo4j connection.
        :type conn: Neo4jConnection
        :raises TypeError: If the connection is not a Neo4jConnection.
        """
        if not isinstance(conn, Neo4jConnection):
            raise TypeError(f"Expected a Neo4jConnection but received {type(conn)}")

        # Set constraints for BioactivityLabel nodes.
        conn.query(
            "CREATE CONSTRAINT FOR (b:BioactivityLabel) REQUIRE (b.name, b.source) IS NODE KEY"
        )
        conn.query("CREATE CONSTRAINT FOR (b:BioactivityLabel) REQUIRE b.name IS NOT NULL")
        conn.query("CREATE CONSTRAINT FOR (b:BioactivityLabel) REQUIRE b.source IS NOT NULL")

    @classmethod
    def create(
        cls,
        conn: Neo4jConnection,
        name: str,
        source: str,
    ) -> None:
        """Create the BioactivityLabel node.

        :param conn: The Neo4j connection.
        :type conn: Neo4jConnection
        :param name: The name of the BioactivityLabel node.
        :type name: str
        :param source: The source of the BioactivityLabel node.
        :type source: str
        :raises TypeError: If the connection is not a Neo4jConnection.
        """
        if not isinstance(conn, Neo4jConnection):
            raise TypeError(f"Expected a Neo4jConnection but received {type(conn)}")

        # Check if any required fields are missing.
        if any([not name, not source]):
            return

        # Check if the BioactivityLabel node already exists.
        if conn.query(
            "MATCH (b:BioactivityLabel {name: $name, source: $source}) RETURN b",
            {"name": name, "source": source},
        ):
            return

        # Create the BioactivityLabel node if it does not exist.
        conn.query(
            "MERGE (b:BioactivityLabel {name: $name, source: $source})",
            {"name": name, "source": source},
        )


class BiosyntheticGeneCluster(Node):
    """The BiosyntheticGeneCluster node."""

    @classmethod
    def set_constraints(cls, conn: Neo4jConnection) -> None:
        """Set constraints on the Neo4j database.

        :param conn: The Neo4j connection.
        :type conn: Neo4jConnection
        :raises TypeError: If the connection is not a Neo4jConnection.
        """
        if not isinstance(conn, Neo4jConnection):
            raise TypeError(f"Expected a Neo4jConnection but received {type(conn)}")

        # Set constraints for BiosyntheticGeneCluster nodes.
        conn.query(
            "CREATE CONSTRAINT FOR (b:BiosyntheticGeneCluster) REQUIRE (b.identifier, b.source) IS NODE KEY"
        )
        conn.query(
            "CREATE CONSTRAINT FOR (b:BiosyntheticGeneCluster) REQUIRE b.identifier IS NOT NULL;"
        )
        conn.query("CREATE CONSTRAINT FOR (b:BiosyntheticGeneCluster) REQUIRE b.source IS NOT NULL")

    @classmethod
    def create(
        cls,
        conn: Neo4jConnection,
        identifier: str,
        source: str,
    ) -> None:
        """Create the BiosyntheticGeneCluster node.

        :param conn: The Neo4j connection.
        :type conn: Neo4jConnection
        :param identifier: The identifier of the BiosyntheticGeneCluster node.
        :type identifier: str
        :param source: The source of the BiosyntheticGeneCluster node.
        :type source: str
        :raises TypeError: If the connection is not a Neo4jConnection.
        """
        if not isinstance(conn, Neo4jConnection):
            raise TypeError(f"Expected a Neo4jConnection but received {type(conn)}")

        # Check if any required fields are missing.
        if any([not identifier, not source]):
            return

        # Check if the BiosyntheticGeneCluster node already exists.
        if conn.query(
            "MATCH (b:BiosyntheticGeneCluster {identifier: $identifier, source: $source}) RETURN b",
            {"identifier": identifier, "source": source},
        ):
            return

        # Create the BiosyntheticGeneCluster node if it does not exist.
        conn.query(
            "MERGE (b:BiosyntheticGeneCluster {identifier: $identifier, source: $source})",
            {"identifier": identifier, "source": source},
        )


class Compound(Node):
    """The Compound node."""

    @classmethod
    def set_constraints(cls, conn: Neo4jConnection) -> None:
        """Set constraints on the Neo4j database.

        :param conn: The Neo4j connection.
        :type conn: Neo4jConnection
        :raises TypeError: If the connection is not a Neo4jConnection.
        """
        if not isinstance(conn, Neo4jConnection):
            raise TypeError(f"Expected a Neo4jConnection but received {type(conn)}")

        # Set constraints for Compound nodes.
        conn.query(
            "CREATE CONSTRAINT FOR (c:Compound) REQUIRE (c.identifier, c.source) IS NODE KEY"
        )
        conn.query("CREATE CONSTRAINT FOR (c:Compound) REQUIRE c.identifier IS NOT NULL")
        conn.query("CREATE CONSTRAINT FOR (c:Compound) REQUIRE c.source IS NOT NULL")
        conn.query("CREATE CONSTRAINT FOR (c:Compound) REQUIRE c.inchikey_connectivity IS NOT NULL")
        conn.query("CREATE CONSTRAINT FOR (c:Compound) REQUIRE c.inchikey IS NOT NULL")
        conn.query("CREATE CONSTRAINT FOR (c:Compound) REQUIRE c.inchi IS NOT NULL")

    @classmethod
    def create(
        cls,
        conn: Neo4jConnection,
        identifier: str,
        source: str,
        inchikey_connectivity: str,
        inchikey: str,
        inchi: str,
    ) -> None:
        """Create the Compound node.

        :param conn: The Neo4j connection.
        :type conn: Neo4jConnection
        :param identifier: The identifier of the Compound node.
        :type identifier: str
        :param source: The source of the Compound node.
        :type source: str
        :param inchikey_connectivity: The InChIKey connectivity of the Compound node.
        :type inchikey_connectivity: str
        :param inchikey: The InChIKey of the Compound node.
        :type inchikey: str
        :param inchi: The InChI of the Compound node.
        :type inchi: str
        :raises TypeError: If the connection is not a Neo4jConnection.
        """
        if not isinstance(conn, Neo4jConnection):
            raise TypeError(f"Expected a Neo4jConnection but received {type(conn)}")

        # Check if any required fields are missing.
        if any([not identifier, not source, not inchikey_connectivity, not inchikey, not inchi]):
            return

        # Check if the Compound node already exists.
        if conn.query(
            "MATCH (c:Compound {identifier: $identifier, source: $source}) RETURN c",
            {"identifier": identifier, "source": source},
        ):
            return

        # Create the Compound node.
        conn.query(
            (
                "MERGE (c:Compound {"
                "identifier: $identifier, "
                "source: $source, "
                "inchikey_connectivity: $inchikey_connectivity, "
                "inchikey: $inchikey, "
                "inchi: $inchi "
                "})"
            ),
            {
                "identifier": identifier,
                "source": source,
                "inchikey_connectivity": inchikey_connectivity,
                "inchikey": inchikey,
                "inchi": inchi,
            },
        )


class Motif(Node):
    """The Motif node."""

    @classmethod
    def set_constraints(cls, conn: Neo4jConnection) -> None:
        """Set constraints on the Neo4j database.

        :param conn: The Neo4j connection.
        :type conn: Neo4jConnection
        :raises TypeError: If the connection is not a Neo4jConnection.
        """
        if not isinstance(conn, Neo4jConnection):
            raise TypeError(f"Expected a Neo4jConnection but received {type(conn)}")

        # Set constraints for Motif nodes.
        conn.query("CREATE CONSTRAINT FOR (m:Motif) REQUIRE m.calculated IS NOT NULL")
        conn.query("CREATE CONSTRAINT FOR (m:Motif) REQUIRE m.calculated IS :: BOOLEAN")
        conn.query(
            "CREATE CONSTRAINT FOR (m:Motif) REQUIRE m.motif_type IS IN ['polyketide', 'peptide']"
        )
        # conn.query("CREATE CONSTRAINT FOR (m:Motif) REQUIRE m.polyketide_type IS NOT NULL IF m.motif_type = 'polyketide'")
        # conn.query("CREATE CONSTRAINT FOR (m:Motif) REQUIRE m.polyketide_decoration_type IS NOT NULL IF m.motif_type = 'polyketide'")
        # conn.query("CREATE CONSTRAINT FOR (m:Motif) REQUIRE m.polyketide_decoration_type IS :: INT IF m.motif_type = 'polyketide'")
        # conn.query("CREATE CONSTRAINT FOR (m:Motif) REQUIRE m.peptide_source IS NOT NULL IF m.motif_type = 'peptide'")
        # conn.query("CREATE CONSTRAINT FOR (m:Motif) REQUIRE m.peptide_cid IS NOT NULL IF m.motif_type = 'peptide'")

    @classmethod
    def create_sub_query(
        cls,
        index: int,
        motif_type: str,
        calculated: bool,
        properties: ty.Dict[str, ty.Any],
    ) -> ty.Tuple[str, ty.Dict[str, ty.Any]]:
        """Create the Motif node.

        :param index: The index of the motif in the motif code.
        :type index: int
        :param motif_type: The type of the Motif node.
        :type motif_type: str
        :param calculated: Whether the Motif node is calculated.
        :type calculated: bool
        :param properties: The properties of the Motif node.
        :type properties: ty.Dict[str, ty.Any]
        :return: The Motif node sub-query and parameters.
        :rtype: ty.Tuple[str, ty.Dict[str, ty.Any]]
        :raises TypeError: If the connection is not a Neo4jConnection.
        :raises ValueError: If any required fields are missing.
        :raises AssertionError: If the motif_type is not one of 'polyketide' or 'peptide'.
        :raises ValueError: If the motif_type is 'polyketide' and the properties are missing
            'polyketide_type' or 'polyketide_decoration_type'.
        :raises ValueError: If the motif_type is 'peptide' and the properties are missing
            'peptide_source' or 'peptide_cid'.
        """
        # If the motif_type is 'polyketide' properties must contain 'polyketide_type' and 'polyketide_decoration_type'.
        if motif_type == "polyketide":

            polyketide_decoration_type = properties.get("polyketide_decoration_type", None)
            if polyketide_decoration_type is None or polyketide_decoration_type == "Any":
                polyketide_decoration_type_query_str = f"polyketide_decoration_type: 'Any'"
            else:
                polyketide_decoration_type_query_str = f"polyketide_decoration_type: {polyketide_decoration_type}"
            
            # Create the polyketide Motif node.
            return (
                f"(u{index + 1}:Motif {{"
                f"motif_type: '{motif_type}', "
                f"calculated: {calculated}, "
                f"polyketide_type: '{properties.get('polyketide_type', None)}', "
                f"{polyketide_decoration_type_query_str}"
                f"}})"
            )

        # If the motif_type is 'peptide' properties must contain 'peptide_source' and 'peptide_cid'.
        elif motif_type == "peptide":
            # Create the peptide Motif node.
            return (
                f"(u{index + 1}:Motif {{"
                f"motif_type: '{motif_type}', "
                f"calculated: {calculated}, "
                f"peptide_source: '{properties.get('peptide_source', None)}', "
                f"peptide_cid: '{properties.get('peptide_cid', None)}'"
                f"}})"
            )
        else:
            raise ValueError(f"Unknown motif type: {motif_type}")

    @classmethod
    def from_string(cls, index: int, motif: str) -> str:
        """Parse motif string (e.g., 'polyketide|A1' or 'peptide|pubchem|1234').

        :param index: The index of the motif in the motif code.
        :type index: int
        :param motif: The motif string.
        :type motif: str
        :return: The Motif node sub-query and parameters.
        :rtype: str
        """
        if match := re.match(r"polyketide\|([A-D])(\d{1,2})", motif):
            motif_type = "polyketide"
            calculated = True
            polyketide_type = match.group(1)
            polyketide_decoration_type = int(match.group(2))

            return Motif.create_sub_query(
                index=index,
                motif_type=motif_type,
                calculated=calculated,
                properties={
                    "polyketide_type": polyketide_type,
                    "polyketide_decoration_type": polyketide_decoration_type,
                },
            )

        elif match := re.match(r"peptide\|(\w+)\|(.+)", motif):
            motif_type = "peptide"
            calculated = True
            peptide_source = match.group(1)
            peptide_cid = match.group(2)

            return Motif.create_sub_query(
                index=index,
                motif_type=motif_type,
                calculated=calculated,
                properties={
                    "peptide_source": peptide_source,
                    "peptide_cid": peptide_cid,
                },
            )

        else:
            raise ValueError(f"Unknown motif: {motif}")


class MotifCode(Node):
    """The MotifCode node."""

    @classmethod
    def set_constraints(cls, conn: Neo4jConnection) -> None:
        """Set constraints on the Neo4j database.

        :param conn: The Neo4j connection.
        :type conn: Neo4jConnection
        :raises TypeError: If the connection is not a Neo4jConnection.
        """
        if not isinstance(conn, Neo4jConnection):
            raise TypeError(f"Expected a Neo4jConnection but received {type(conn)}")

        # Set constraints for MotifCode nodes.
        conn.query(
            "CREATE CONSTRAINT FOR (m:MotifCode) REQUIRE (m.compound_identifier, m.compound_source, m.calculated) IS NODE KEY"
        )  # noqa: E501
        conn.query(
            "CREATE CONSTRAINT FOR (m:MotifCode) REQUIRE m.compound_identifier IS NOT NULL"
        )  # noqa: E501
        conn.query(
            "CREATE CONSTRAINT FOR (m:MotifCode) REQUIRE m.compound_source IS NOT NULL"
        )  # noqa: E501
        conn.query(
            "CREATE CONSTRAINT FOR (m:MotifCode) REQUIRE m.calculated IS NOT NULL"
        )  # noqa: E501
        conn.query(
            "CREATE CONSTRAINT FOR (m:MotifCode) REQUIRE m.calculated IS :: BOOLEAN"
        )  # noqa: E501
        # conn.query("CREATE CONSTRAINT FOR (m:MotifCode) REQUIRE (m.src IS NOT NULL) IF m.calculated")  # noqa: E501
        # conn.query("CREATE CONSTRAINT FOR (m:MotifCode) REQUIRE (m.src IS NULL) IF NOT m.calculated")  # noqa: E501
        # conn.query("CREATE CONSTRAINT FOR (m:MotifCode) REQUIRE (m.applied_reactions IS NOT NULL) IF m.calculated")  # noqa: E501
        # conn.query("CREATE CONSTRAINT FOR (m:MotifCode) REQUIRE (m.applied_reactions IS NULL) IF NOT m.calculated")  # noqa: E501
        # conn.query("CREATE CONSTRAINT FOR (m:MotifCode) REQUIRE all(x IN m.applied_reactions WHERE x IS NOT NULL) IF m.calculated")  # noqa: E501

    @classmethod
    def create_sub_query(
        cls,
        compound_identifier: str,
        compound_source: str,
        calculated: bool,
        applied_reactions: ty.Optional[ty.List[str]] = None,
        src: ty.Optional[str] = None,
    ) -> ty.Tuple[str, ty.Dict[str, ty.Any]]:
        """Create the MotifCode node.

        :param compound_identifier: The identifier of the MotifCode node.
        :type compound_identifier: str
        :param compound_source: The source of the MotifCode node.
        :type compound_source: str
        :param calculated: Whether the MotifCode node is calculated.
        :type calculated: bool
        :param applied_reactions: The applied reactions of the MotifCode node.
        :type applied_reactions: ty.Optional[ty.List[str]]
        :param src: The source of the MotifCode node.
        :type src: str
        :return: The MotifCode node sub-query and parameters.
        :rtype: ty.Tuple[str, ty.Dict[str, ty.Any]]
        :raises TypeError: If the connection is not a Neo4jConnection.
        """
        return (
            (
                "(m:MotifCode {"
                "compound_identifier: $compound_identifier, "
                "compound_source: $compound_source, "
                "calculated: $calculated, "
                "applied_reactions: $applied_reactions, "
                "src: $src"
                "})"
            ),
            {
                "compound_identifier": compound_identifier,
                "compound_source": compound_source,
                "calculated": calculated,
                "applied_reactions": applied_reactions,
                "src": src,
            },
        )


class Organism(Node):
    """The Organism node."""

    @classmethod
    def set_constraints(cls, conn: Neo4jConnection) -> None:
        """Set constraints on the Neo4j database.

        :param conn: The Neo4j connection.
        :type conn: Neo4jConnection
        :raises TypeError: If the connection is not a Neo4jConnection.
        """
        if not isinstance(conn, Neo4jConnection):
            raise TypeError(f"Expected a Neo4jConnection but received {type(conn)}")

        # Set constraints for Organism nodes.
        conn.query("CREATE CONSTRAINT FOR (o:Organism) REQUIRE o.ncbi_id IS NODE KEY")
        conn.query("CREATE CONSTRAINT FOR (o:Organism) REQUIRE o.ncbi_id IS NOT NULL")
        conn.query("CREATE CONSTRAINT FOR (o:Organism) REQUIRE o.type IS NOT NULL")
        conn.query("CREATE CONSTRAINT FOR (o:Organism) REQUIRE o.genus IS NOT NULL")
        conn.query("CREATE CONSTRAINT FOR (o:Organism) REQUIRE o.species IS NOT NULL")

    @classmethod
    def create(
        cls,
        conn: Neo4jConnection,
        ncbi_id: str,
        type: str,
        genus: str,
        species: str,
    ) -> None:
        """Create the Organism node.

        :param conn: The Neo4j connection.
        :type conn: Neo4jConnection
        :param ncbi_id: The NCBI ID of the Organism node.
        :type ncbi_id: str
        :param type: The type of the Organism node.
        :type type: str
        :param genus: The genus of the Organism node.
        :type genus: str
        :param species: The species of the Organism node.
        :type species: str
        :raises TypeError: If the connection is not a Neo4jConnection.
        """
        if not isinstance(conn, Neo4jConnection):
            raise TypeError(f"Expected a Neo4jConnection but received {type(conn)}")

        # Check if any required fields are missing.
        if any([not ncbi_id, not type, not genus, not species]):
            return

        # Check if the Organism node already exists.
        if conn.query("MATCH (o:Organism {ncbi_id: $ncbi_id}) RETURN o", {"ncbi_id": ncbi_id}):
            return

        # Create the Organism node if it does not exist.
        conn.query(
            (
                "MERGE (o:Organism {"
                "ncbi_id: $ncbi_id, "
                "type: $type, "
                "genus: $genus, "
                "species: $species "
                "})"
            ),
            {
                "ncbi_id": ncbi_id,
                "type": type,
                "genus": genus,
                "species": species,
            },
        )


class Pathway(Node):
    """The Pathway node."""

    @classmethod
    def set_constraints(cls, conn: Neo4jConnection) -> None:
        """Set constraints on the Neo4j database.

        :param conn: The Neo4j connection.
        :type conn: Neo4jConnection
        :raises TypeError: If the connection is not a Neo4jConnection.
        """
        if not isinstance(conn, Neo4jConnection):
            raise TypeError(f"Expected a Neo4jConnection but received {type(conn)}")

        # Set constraints for Pathway nodes.
        conn.query("CREATE CONSTRAINT FOR (p:Pathway) REQUIRE (p.name, p.source) IS NODE KEY")
        conn.query("CREATE CONSTRAINT FOR (p:Pathway) REQUIRE p.name IS NOT NULL")
        conn.query("CREATE CONSTRAINT FOR (p:Pathway) REQUIRE p.source IS NOT NULL")

    @classmethod
    def create(
        cls,
        conn: Neo4jConnection,
        name: str,
        source: str,
    ) -> None:
        """Create the Pathway node.

        :param conn: The Neo4j connection.
        :type conn: Neo4jConnection
        :param name: The name of the Pathway node.
        :type name: str
        :param source: The source of the Pathway node.
        :type source: str
        :raises TypeError: If the connection is not a Neo4jConnection.
        """
        if not isinstance(conn, Neo4jConnection):
            raise TypeError(f"Expected a Neo4jConnection but received {type(conn)}")

        # Check if any required fields are missing.
        if any([not name, not source]):
            return

        # Check if the Pathway node already exists.
        if conn.query(
            "MATCH (p:Pathway {name: $name, source: $source}) RETURN p",
            {"name": name, "source": source},
        ):
            return

        # Create the Pathway node if it does not exist.
        conn.query(
            "MERGE (p:Pathway {name: $name, source: $source})", {"name": name, "source": source}
        )


class Protocluster(Node):
    """The Protocluster node."""

    @classmethod
    def set_constraints(cls, conn: Neo4jConnection) -> None:
        """Set constraints on the Neo4j database.

        :param conn: The Neo4j connection.
        :type conn: Neo4jConnection
        :raises TypeError: If the connection is not a Neo4jConnection.
        """
        if not isinstance(conn, Neo4jConnection):
            raise TypeError(f"Expected a Neo4jConnection but received {type(conn)}")

        # Set constraints for Protocluster nodes.
        conn.query(
            "CREATE CONSTRAINT FOR (p:Protocluster) REQUIRE (p.identifier, p.source) IS NODE KEY"
        )
        conn.query("CREATE CONSTRAINT FOR (p:Protocluster) REQUIRE p.identifier IS NOT NULL")
        conn.query("CREATE CONSTRAINT FOR (p:Protocluster) REQUIRE p.source IS NOT NULL")

    @classmethod
    def create(
        cls,
        conn: Neo4jConnection,
        identifier: str,
        source: str,
    ) -> None:
        """Create the Protocluster node.

        :param conn: The Neo4j connection.
        :type conn: Neo4jConnection
        :param identifier: The identifier of the Protocluster node.
        :type identifier: str
        :param source: The source of the Protocluster node.
        :type source: str
        :raises TypeError: If the connection is not a Neo4jConnection.
        """
        if not isinstance(conn, Neo4jConnection):
            raise TypeError(f"Expected a Neo4jConnection but received {type(conn)}")

        # Check if any required fields are missing.
        if any([not identifier, not source]):
            return

        # Create the Protocluster node if it does not exist.
        conn.query(
            "MERGE (p:Protocluster {identifier: $identifier, source: $source})",
            {"identifier": identifier, "source": source},
        )
