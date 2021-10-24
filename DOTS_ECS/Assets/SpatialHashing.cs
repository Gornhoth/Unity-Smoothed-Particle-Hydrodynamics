using UnityEngine;

public static class SpatialHashing
{
    public static float CellSize = 1f;
    public static int Dimensions = 10;

    public static Vector3Int GetCell(Vector3 position, float cellSize)
    {
        return new Vector3Int((int)( position.x / cellSize ), (int)( position.y / cellSize ), (int)( position.z / cellSize ));
    }

    public static int Hash(Vector3Int cell, int dimensions)
    {
        return cell.x + dimensions * ( cell.y + dimensions * cell.z );
    }
}