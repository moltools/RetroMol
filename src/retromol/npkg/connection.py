# -*- coding: utf-8 -*-

"""This module contains the connector for interacting with the Neo4j database."""

import logging
import typing as ty

from neo4j import GraphDatabase


class Neo4jConnection:
    """The connector for interacting with the Neo4j database.

    Adapted from: https://towardsdatascience.com/create-a-graph-database-in-neo4j-using-python-4172d40f89c4
    """

    def __init__(self, uri: str, user: str, pwd: str) -> None:
        """Initialize the Neo4j connection.

        :param uri: The URI of the Neo4j database.
        :type uri: str
        :param user: The user of the Neo4j database.
        :type user: str
        :param pwd: The password of the Neo4j database.
        :type pwd: str
        :raises Exception: If the driver cannot be created.
        """
        self.logger = logging.getLogger(__name__)
        self._uri = uri
        self._user = user
        self._pwd = pwd
        self._driver = None

        try:
            self._driver = GraphDatabase.driver(
                uri=self._uri,
                auth=(self._user, self._pwd),
            )
        except Exception as e:
            msg = f"Failed to create the driver: {e}"
            self.logger.error(msg)
            raise e

    def close(self) -> None:
        """Close the Neo4j connection."""
        if self._driver is not None:
            self._driver.close()

    def query(
        self,
        query: str,
        parameters: ty.Dict[str, ty.Any] = None,
        db: str = None,
        batch_size: int = None,
    ) -> ty.List[ty.Any]:
        """Execute a query on the Neo4j database.

        :param query: The query to execute.
        :type query: str
        :param parameters: The parameters of the query.
        :type parameters: ty.Dict[str, ty.Any]
        :param db: The database to query.
        :type db: str
        :param batch_size: The batch size of the query.
        :type batch_size: int
        :return: The result of the query.
        :rtype: ty.Any
        :raises Exception: If the query fails.
        """
        assert self._driver is not None, "Driver not initialized!"

        session = None
        response = None

        try:
            session = self._driver.session(database=db) if db else self._driver.session()
            response = list(session.run(query, parameters, batch_size=batch_size))

        except Exception as e:
            msg = f"Failed to execute query: {e}"
            self.logger.error(msg)
            raise e

        finally:
            if session is not None:
                session.close()

        return response
